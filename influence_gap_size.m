% Calculation of an interval power spectrum from signals with missing data
%
% Investigation of the influence of the gap size of the missing data
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

n_md = 1; % number of missing data
l_md = 1; % length of the missing data gap
dist_md = 'uniform'; % distribution of the missing data in the signal

count = 1; % loop-internal counter

gap = 0:0.5:31.5;
S_signal_std = cell(1,length(gap));
S_signal_minmax = cell(1,length(gap));

S_std = cell(1,length(gap));
S_minmax = cell(1,length(gap));

S_area_std = zeros(1,length(gap));
S_area_minmax= zeros(1,length(gap));

%% generate missing data
for gap_it = gap

    % get missing data values and position (start in the middle of the
    % signal and let the gap "grow" in each for-loop)
    missing_data = signal(round(n_signal/2-gap_it):round(n_signal/2+gap_it));
    md_pos = round(n_signal/2-gap_it):round(n_signal/2+gap_it);

    % rearrange missing data gaps
    missing_data_coord = [reshape(md_pos',1,numel(md_pos)); reshape(missing_data',1,numel(missing_data))];
    [md, sort_ind] = sort(missing_data_coord(1,:));
    missing_data_coord(2,:) = missing_data_coord(2,sort_ind);

    % reconstruct missing data by adding interval uncertainty x0 to true
    % data point
    int_signal_reconstructed_std = [signal;signal];
    int_signal_reconstructed_std(:,md) = [int_signal_reconstructed_std(1,md)-x0; int_signal_reconstructed_std(2,md)+x0];

    % reconstruct missing data by intervals with signal widths
    s_max = max(signal);
    s_min = min(signal);
    int_signal_reconstructed_signalminmax = [signal;signal];
    int_signal_reconstructed_signalminmax(:,md) = ones(size(md)).*[s_min; s_max];

    S_signal_std{count} = int_signal_reconstructed_std;
    S_signal_minmax{count} = int_signal_reconstructed_signalminmax;

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
    S_rec_interval_minmax = periodogram_without_fft([S_rec_interval_minmax(:,1) S_rec_interval_minmax]);

    S_std{count} = S_rec_interval_std;
    S_minmax{count} = S_rec_interval_minmax;

    S_area_std(count) = sum(S_rec_interval_std(2,:) - S_rec_interval_std(1,:));
    S_area_minmax(count) = sum(S_rec_interval_minmax(2,:) - S_rec_interval_minmax(1,:));

    count = count + 1;
end

%% plot area between bounds

kk_new = 1:count-1;
figure; hold on; grid on;
plot(kk_new, smooth(S_area_std))
plot(kk_new, smooth(S_area_minmax))
xlabel('Gap size')
ylabel('Area between bounds')
legend('Reconstruction method (1)','Reconstruction method (2)', 'Location', 'northwest')
xlim([0 64])

%% plots reconstruction method (1)

figure;
subplot(4,2,2); hold on;
plot_intervalPSD(w, S_std{1});
plot(w, S_target(1:Nw), 'Color', [0.1216 0.4667 0.7059], 'LineWidth', 1);
xlim([0 w(end)])
xlabel('Frequency (rad/s)')
ylabel('Spectral densitiy')

subplot(4,2,1); hold on;
plot_intervalsignal(t, S_signal_std{1})
xlim([0 t(end)])
xlabel('Time (s)')
ylabel('Wave height (m)')

subplot(4,2,4); hold on;
plot_intervalPSD(w, S_std{20});
plot(w, S_target(1:Nw), 'Color', [0.1216 0.4667 0.7059], 'LineWidth', 1);
xlim([0 w(end)])
xlabel('Frequency (rad/s)')
ylabel('Spectral densitiy')

subplot(4,2,3); hold on;
plot_intervalsignal(t, S_signal_std{20})
xlim([0 t(end)])
xlabel('Time (s)')
ylabel('Wave height (m)')

subplot(4,2,6); hold on;
plot_intervalPSD(w, S_std{40});
plot(w, S_target(1:Nw), 'Color', [0.1216 0.4667 0.7059], 'LineWidth', 1);
xlim([0 w(end)])
xlabel('Frequency (rad/s)')
ylabel('Spectral densitiy')

subplot(4,2,5); hold on;
plot_intervalsignal(t, S_signal_std{40})
xlim([0 t(end)])
xlabel('Time (s)')
ylabel('Wave height (m)')

subplot(4,2,8); hold on;
plot_intervalPSD(w, S_std{60});
plot(w, S_target(1:Nw), 'Color', [0.1216 0.4667 0.7059], 'LineWidth', 1);
xlim([0 w(end)])
xlabel('Frequency (rad/s)')
ylabel('Spectral densitiy')

subplot(4,2,7); hold on;
plot_intervalsignal(t, S_signal_std{60})
xlim([0 t(end)])
xlabel('Time (s)')
ylabel('Wave height (m)')

%% plots reconstruction method (2)

figure;
subplot(4,2,2); hold on;
plot_intervalPSD(w, S_minmax{1});
plot(w, S_target(1:Nw), 'Color', [0.1216 0.4667 0.7059], 'LineWidth', 1);
xlim([0 w(end)])
xlabel('Frequency (rad/s)')
ylabel('Spectral densitiy')

subplot(4,2,1); hold on;
plot_intervalsignal(t, S_signal_minmax{1})
xlim([0 t(end)])
xlabel('Time (s)')
ylabel('Wave height (m)')

subplot(4,2,4); hold on;
plot_intervalPSD(w, S_minmax{20});
plot(w, S_target(1:Nw), 'Color', [0.1216 0.4667 0.7059], 'LineWidth', 1);
xlim([0 w(end)])
xlabel('Frequency (rad/s)')
ylabel('Spectral densitiy')

subplot(4,2,3); hold on;
plot_intervalsignal(t, S_signal_minmax{20})
xlim([0 t(end)])
xlabel('Time (s)')
ylabel('Wave height (m)')

subplot(4,2,6); hold on;
plot_intervalPSD(w, S_minmax{40});
plot(w, S_target(1:Nw), 'Color', [0.1216 0.4667 0.7059], 'LineWidth', 1);
xlim([0 w(end)])
xlabel('Frequency (rad/s)')
ylabel('Spectral densitiy')

subplot(4,2,5); hold on;
plot_intervalsignal(t, S_signal_minmax{40})
xlim([0 t(end)])
xlabel('Time (s)')
ylabel('Wave height (m)')

subplot(4,2,8); hold on;
plot_intervalPSD(w, S_minmax{60});
plot(w, S_target(1:Nw), 'Color', [0.1216 0.4667 0.7059], 'LineWidth', 1);
xlim([0 w(end)])
xlabel('Frequency (rad/s)')
ylabel('Spectral densitiy')

subplot(4,2,7); hold on;
plot_intervalsignal(t, S_signal_minmax{60})
xlim([0 t(end)])
xlabel('Time (s)')
ylabel('Wave height (m)')
