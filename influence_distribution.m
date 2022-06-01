% Calculation of an interval power spectrum from signals with missing data
%
% Investigation of the influence of the distribution of the missing data 
% within the signal
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
n_signal = length(signal);
x0 = std(signal); % interval width is sample standard deviation of the signal

periodogram_without_fft = @(x) abs(x).^2.*dt^2/T/(2*pi); % periodogram to estimate the PSD of a signal

count = 1; % loop-internal counter

l_md = 1; % length of the missing data gap
n_md = [4 8 16 32]; % number of missing data gaps
S_std = cell(1,2*length(n_md));

%% generate missing data
for n_md_it = n_md
    for dist_md = ["uniform" "binomial"]

        % get missing data values and position (each row in the variable
        % "missing_data" and "md_pos" corresponds to one missing data gap)
        [missing_data, md_pos] = generate_missing_data(signal, n_md_it, l_md, dist_md);
        
        % rearrange missing data gaps for internal calculations
        missing_data_coord = [reshape(md_pos',1,numel(md_pos)); reshape(missing_data',1,numel(missing_data))];
        [md, sort_ind] = sort(missing_data_coord(1,:));
        missing_data_coord(2,:) = missing_data_coord(2,sort_ind);

        % reconstruct missing data by adding interval uncertainty x0 to
        % true data point
        int_signal_reconstructed_std = [signal;signal];
        int_signal_reconstructed_std(:,md) = [int_signal_reconstructed_std(1,md)-x0; int_signal_reconstructed_std(2,md)+x0];

        % reconstruct missing data by intervals with signal widths
        s_max = max(signal);
        s_min = min(signal);
        int_signal_reconstructed_signalminmax = [signal;signal];
        int_signal_reconstructed_signalminmax(:,md) = ones(size(md)).*[s_min; s_max];
        
        %% compute interval PSDs
        
        % target spectrum
        S_target = 0;
        k = 0:Nw-1;
        for n=1:Nt
            S_target = S_target + exp(-1i*2*pi*(n-1)*k/Nt)*signal(n);
        end
        S_target = periodogram_without_fft(S_target);

        % interval spectrum of signal_x0
        S_rec_interval_std = zeros(2,Nt/2-1);
        S_rec_interval_minmax = zeros(2,Nt/2-1);
        for freq = 1:Nt/2-1
            S_rec_interval_std(:,freq) = interval_DFT(int_signal_reconstructed_std,freq);
            S_rec_interval_minmax(:,freq) = interval_DFT(int_signal_reconstructed_signalminmax,freq);
        end
        % for numerical issues it is not possible to calculate the convex hull of
        % the first frequency k. Therefore the first frequency is added
        % artificially
        S_rec_interval_std = periodogram_without_fft([S_rec_interval_std(:,1) S_rec_interval_std]);

        S_std{count} = S_rec_interval_std;

        count = count + 1;
    end
end

%% plots
figure; subplot(4,2,1); hold on;
plot_intervalPSD(w, S_std{1});
plot(w, S_target(1:Nw), 'Color', [0.1216 0.4667 0.7059], 'LineWidth', 1);
xlim([0 w(end)])
xlabel('Frequency (rad/s)')
ylabel('Spectral densitiy')
title('Uniformly distributed')

subplot(4,2,2); hold on;
plot_intervalPSD(w, S_std{2});
plot(w, S_target(1:Nw), 'Color', [0.1216 0.4667 0.7059], 'LineWidth', 1);
xlim([0 w(end)])
xlabel('Frequency (rad/s)')
ylabel('Spectral densitiy')
title('Binomially distributed')

subplot(4,2,3); hold on;
plot_intervalPSD(w, S_std{3});
plot(w, S_target(1:Nw), 'Color', [0.1216 0.4667 0.7059], 'LineWidth', 1);
xlim([0 w(end)])
xlabel('Frequency (rad/s)')
ylabel('Spectral densitiy')

subplot(4,2,4); hold on;
plot_intervalPSD(w, S_std{4});
plot(w, S_target(1:Nw), 'Color', [0.1216 0.4667 0.7059], 'LineWidth', 1);
xlim([0 w(end)])
xlabel('Frequency (rad/s)')
ylabel('Spectral densitiy')

subplot(4,2,5); hold on;
plot_intervalPSD(w, S_std{5});
plot(w, S_target(1:Nw), 'Color', [0.1216 0.4667 0.7059], 'LineWidth', 1);
xlim([0 w(end)])
xlabel('Frequency (rad/s)')
ylabel('Spectral densitiy')

subplot(4,2,6); hold on;
plot_intervalPSD(w, S_std{6});
plot(w, S_target(1:Nw), 'Color', [0.1216 0.4667 0.7059], 'LineWidth', 1);
xlim([0 w(end)])
xlabel('Frequency (rad/s)')
ylabel('Spectral densitiy')

subplot(4,2,7); hold on;
plot_intervalPSD(w, S_std{7});
plot(w, S_target(1:Nw), 'Color', [0.1216 0.4667 0.7059], 'LineWidth', 1);
xlim([0 w(end)])
xlabel('Frequency (rad/s)')
ylabel('Spectral densitiy')

subplot(4,2,8); hold on;
plot_intervalPSD(w, S_std{8});
plot(w, S_target(1:Nw), 'Color', [0.1216 0.4667 0.7059], 'LineWidth', 1);
xlim([0 w(end)])
xlabel('Frequency (rad/s)')
ylabel('Spectral densitiy')



