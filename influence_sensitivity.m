% Calculation of an interval power spectrum from signals with missing data
%
% Investigation of the sensitivity of the reconstructed interval width
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

periodogram_without_fft = @(x) abs(x).^2.*dt^2/T/(2*pi); % periodogram to estimate the PSD of a signal

n_md = 1; % number of missing data
l_md = 1; % length of the missing data gap
dist_md = 'uniform'; % distribution of the missing data in the signal

interval_height = 0.1:0.1:10;
gap = 0:1:5;

S_width = zeros(length(gap),length(interval_height));
S_area = zeros(length(gap),length(interval_height));

%% generate missing data
for gap_it = gap
    count = 1;
    for interval_height_it = interval_height
        
        % get missing data values and position (start in the middle of the 
        % signal and let the gap "grow" in each for-loop)
        missing_data = signal(n_signal/2-gap_it:n_signal/2+gap_it);
        md_pos = n_signal/2-gap_it:n_signal/2+gap_it;

        % rearrange missing data gaps
        missing_data_coord = [reshape(md_pos',1,numel(md_pos)); reshape(missing_data',1,numel(missing_data))];
        [md, sort_ind] = sort(missing_data_coord(1,:));

        % reconstruct missing data by adding interval uncertainty x0 to
        % true data point
        int_signal_reconstructed_std = [signal;signal];
        int_signal_reconstructed_std(:,md) = [int_signal_reconstructed_std(1,md)-interval_height_it; int_signal_reconstructed_std(2,md)+interval_height_it];

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

        S_width(gap_it+1,count) = abs(S_rec_interval_std(1,11)-S_rec_interval_std(2,11));
        S_area(gap_it+1,count) = sum(S_rec_interval_std(2,:) - S_rec_interval_std(1,:));

        count = count + 1;

    end
end

%% plots 
figure; hold on; grid on;
for i = 1:6
    plot(interval_height,S_width(i,:))
end
xlabel('Interval uncertainty \xi')
ylabel('Interval width at \omega_p')
legend('Gap size = 1','Gap size = 3','Gap size = 5','Gap size = 7','Gap size = 9','Gap size = 11','Location','northwest')

figure; hold on; grid on;
for i = 1:6
    plot(interval_height,S_area(i,:))
end
xlabel('Interval uncertainty \xi')
ylabel('Area between bounds')
legend('Gap size = 1','Gap size = 3','Gap size = 5','Gap size = 7','Gap size = 9','Gap size = 11','Location','northwest')

