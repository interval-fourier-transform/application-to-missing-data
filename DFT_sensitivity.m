% Influence of a simple linear interpolation on the amplitude of the DFT.
% The linear reconstruction is obtained by the two adjacent non-missing 
% data points.
%
% Marco Behrendt
% Institute for Risk and Reliability, Leibniz Universität Hannover
% behrendt@irz.uni-hannover.de
% https://github.com/marcobehrendt
%
% Date: 17/03/2022

clear variables; close all; clc;
addpath(genpath(pwd)); 

load('signal_64.mat')
x0 = std(signal);

periodogram = @(x) abs(fft(x)).^2.*dt^2/T/(2*pi); % periodogram to estimate the PSD of a signal

%% Target PSD for comparison
S_target = periodogram(signal);
S_target = S_target(1:length(w));

%% 5% missing dat
pos_md = randperm(Nt);
pos_md = sort(pos_md(1:round(Nt*0.05)));

% generate vector where missing data is represented by NaN
signal1 = signal; 
signal1(pos_md) = NaN;

% reconstruct missing data by linear interpolation of the adjacent points
interp_method = 'linear';
signal1 = fillmissing(signal1, interp_method);

% estimate the PSD
S1 = periodogram(signal1);
S1 = S1(1:length(w));

%% 10% missing dat
pos_md = randperm(Nt);
pos_md = sort(pos_md(1:round(Nt*0.1)));

% generate vector where missing data is represented by NaN
signal2 = signal; 
signal2(pos_md) = NaN;

% reconstruct missing data by linear interpolation
interp_method = 'linear';
signal2 = fillmissing(signal2, interp_method);

S2 = periodogram(signal2);
S2 = S2(1:length(w));

%% 25% missing dat
pos_md = randperm(Nt);
pos_md = sort(pos_md(1:round(Nt*0.25)));

% generate vector where missing data is represented by NaN
signal3 = signal; 
signal3(pos_md) = NaN;

% reconstruct missing data by linear interpolation
interp_method = 'linear';
signal3 = fillmissing(signal3, interp_method);

S3 = periodogram(signal3);
S3 = S3(1:length(w));

%% plot results
figure; hold on; grid on;
plot(w, S_target)
plot(w, S1, '--')
plot(w, S2, ':', 'Color', [0.49 0.18 0.56])
plot(w, S3, '-.', 'Color', [0.47 0.67 0.19])
legend('Target spectrum','Reconstruction with 5% missing data','Reconstruction with 10% missing data','Reconstruction with 25% missing data')
xlabel('Frequency (rad/s)')
ylabel('PSD (m^2s)')