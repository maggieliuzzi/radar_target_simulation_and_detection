clear
clc
close all


%% RADAR SIMULATION


% speed of light: 3e8

%% Radar Specs
% Frequency of operation: 77GHz (quite common in automotive applications)
% Max range: 200m
% Range resolution: 1m
% Max velocity: 100m/s

%% Target Position & Velocity
% initial position/ range
d0 = 80;
% velocity (assumed constant)
v0 = -50;

%% FMCW Waveform Design

d_res = 1;
speed_of_light = 3*10^8;
RMax = 200;

Bsweep = speed_of_light / (2 * d_res);  % bandwidth
Tchirp = 5.5 * 2 * RMax / speed_of_light; % chirp time
alpha = Bsweep / Tchirp; % slope of FMCW chirps
fc = 77e9;              % carrier frequency of radar 
                                                     
Nd = 128;                % number of chirps in one sequence
Nr = 1024;               % number of samples on each chirp    

t = linspace(0, Nd * Tchirp, Nr * Nd); % total time for samples

% vectors for Tx, Rx and Mix based on the total samples input
Tx = zeros(1, length(t));    % transmitted signal
Rx = zeros(1, length(t));    % received signal
Mix = zeros(1, length(t));   % beat signal

% vectors for range covered and time delay
r_t = zeros(1, length(t));
td = zeros(1, length(t));

%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

for i = 1:length(t)         
    % for each timestamp update the range of the target for constant velocity
    r_t(i) = d0 + v0 * t(i);
    td(i) = 2 * r_t(i) / speed_of_light;
    
    % for each time sample update the transmitted and received signal
    Tx(i) = cos(2 * pi * (fc * t(i) + alpha * t(i)^2 / 2));
    Rx(i) = cos(2 * pi * (fc * (t(i) - td(i)) + (alpha * (t(i) - td(i))^2) / 2));
    
    % now by mixing the transmit and receive generate the beat signal by performing
    % element-wise matrix multiplication of transmit and receiver signal
    Mix(i) = Tx(i) .* Rx(i);  % beat signal
    
end


%% RANGE MEASUREMENT


% reshape vector into Nr*Nd array. Nr and Nd here also define the size of range and doppler FFT respectively

% run the Fast Fourier Transform (FFT) on the beat signal along the range bins dimension (Nr) and normalise
sig_fft = fft(Mix, Nr) ./ Nr;

% take absolute value of FFT output
sig_fft = abs(sig_fft);       

% Output of FFT is double-sided signal, but we are interested in one side of the spectrum only, hence we throw out half of the samples
sig_fft = sig_fft(1 : (Nr / 2));

figure ('Name', 'Range from first FFT')  % plot range
plot(sig_fft); grid minor                % plot FFT output 
axis([0 200 0 1]);
xlabel('Measured range');


%% RANGE DOPPLER RESPONSE

% Run a 2D FFT on the mixed signal (beat signal) output and generate a Range Doppler Map (RDM)

% The output of the 2D FFT is an image that has reponse in the range and doppler FFT bins. 
% So, it is important to convert the axis from bin sizes to range and doppler based on their max values

Mix = reshape(Mix, [Nr, Nd]);

% 2D FFT using the FFT size for both dimensions
sig_fft2 = fft2(Mix, Nr, Nd);

% Taking just one side of signal from range dimension.
sig_fft2 = sig_fft2(1 : Nr / 2, 1 : Nd);
sig_fft2 = fftshift(sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM);

% use surf function to plot output of 2D-FFT and show axis in both dimensions
doppler_axis = linspace(-100, 100, Nd);
range_axis = linspace(-200, 200, Nr/2) * ((Nr / 2) / 400);
figure('Name', 'Range Doppler Map'), surf(doppler_axis, range_axis, RDM);


%% CFAR Implementation

% slide window through the complete Range Doppler Map

% selecting number of training cells in both the dimensions
Tcr = 10;
Tcd = 4;

% selecting number of guard cells in both dimensions around the Cell Under Test (CUT) for accurate estimation
Gcr = 5;
Gcd = 2;

% offset the threshold by signal-to-noise-ratio (SNR) value in dB
offset = 1.4;

% create vector to store noise level for each iteration on training cells
noise_level = zeros(Nr / 2 - 2 * (Tcd + Gcd), Nd - 2 * (Tcr + Gcr));
gridSize = (2 * Tcr + 2 * Gcr + 1) * (2 * Tcd + 2 * Gcd + 1);
trainingCellsNum = gridSize - (2 * Gcr + 1) * (2 * Gcd + 1);

% designing loop - slide the CUT across range doppler map by giving margins 
% at the edges for training and guard cells.
% For every iteration sum the signal level within all the training cells
% To sum convert the value from logarithmic to linear using db2pow function
% Average the summed values for all of the training cells used
% After averaging, convert it back to logarithimic using pow2db
% Add the offset to it to determine the threshold 
% Next, compare the signal under CUT with this threshold:
% if the CUT level > threshold, assign it a value of 1, else equate it to 0

CFAR_sig = zeros(size(RDM));

% Use RDM[x,y] as the matrix from the output of 2D FFT for implementing CFAR
for j = 1 : Nd - 2 * (Tcr + Gcr)
    for i = 1 : Nr / 2 - 2 * (Tcd + Gcd)
        % to extract only the training cells, first get the sliding patch and, 
        % after converting it to power from decibel, set to zero whatever isn't 
        % in the position of training cells;
        % the zero submatrix will not contribute to the sum over the whole patch and, 
        % by dividing for the number of training cells I will get the noise mean level
        
        trainingCellsPatch = db2pow(RDM( i : i + 2 * (Tcd + Gcd), j : j + 2 * (Gcr + Tcr)));
        trainingCellsPatch(Tcd + 1 : end - Tcd, Tcr + 1 : end - Tcr) = 0;
        
        noise_level(i,j) = pow2db(sum(sum(trainingCellsPatch))/trainingCellsNum);
        sigThresh = noise_level(i, j) * offset;
        if RDM(i + (Tcd + Gcd), j + (Tcd + Gcr)) > sigThresh
            CFAR_sig(i + (Tcd + Gcd), j + (Tcd + Gcr)) = 1;
        else
            CFAR_sig(i + (Tcd + Gcd), j + (Tcd + Gcr)) = 0;
        end
    end
end

% The process above will generate a thresholded block, which is smaller than the 
% Range Doppler Map as the CUT cannot be located at the edges of matrix
% Hence, few cells will not be thresholded
% To keep the map size same set those values to 0  %% necessary?

% display CFAR output using the surf function
figure('Name', 'CA-CFAR Filtered RDM'), surf(doppler_axis, range_axis, CFAR_sig);
colorbar;
