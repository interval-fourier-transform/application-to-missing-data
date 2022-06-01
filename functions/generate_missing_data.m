function [missing_data, missing_data_pos] = generate_missing_data(signal, n_md, l_md, dist_md)
%Generation of missing data for a given signal
% This function generates randomly missing data in a given signal. Input
% parameters are the signal itself, the number of missing data n_md, the
% length l_md of each missing data gap and the distribution dist_md, i.e.,
% how the missing data gap is distributed within the signal.
%
% Author:
% Marco Behrendt
% Institute for Risk and Reliability, Leibniz Universität Hannover
% behrendt@irz.uni-hannover.de
% https://github.com/marcobehrendt
%
% Date: 11/03/2022

Nt = length(signal);

missing_data = zeros(n_md, l_md);
missing_data_pos = zeros(n_md, l_md);

if rem(l_md,2) == 1 % if l_md is odd
    for k = 1:n_md
        md_pos(k) = new_rnd(Nt,dist_md);
        seq(k,:) = md_pos(k)-floor(l_md/2):md_pos(k)+floor(l_md/2);

        if k == 1
            while any(seq(k,:) > Nt) || any(seq(k,:) < 1)
                md_pos(k) = new_rnd(Nt,dist_md);
                seq(k,:) = md_pos(k)-floor(l_md/2):md_pos(k)+floor(l_md/2);
            end
        else
            while any(intersect(seq(1:k-1, :), seq(k,:))) || any(seq(k,:) > Nt) || any(seq(k,:) < 1)
                md_pos(k) = new_rnd(Nt,dist_md);
                seq(k,:) = md_pos(k)-floor(l_md/2):md_pos(k)+floor(l_md/2);
            end
        end
    end

    for k = 1:n_md
        missing_data(k,:) = signal(seq(k,:));
        missing_data_pos(k,:) = seq(k,:);
    end

else % if l_md is even

    for k = 1:n_md
        md_pos(k) = new_rnd(Nt,dist_md);
        seq(k,:) = md_pos(k)-(l_md/2-1):md_pos(k)+(l_md/2);

        if k == 1
            while any(seq(k,:) > Nt) || any(seq(k,:) < 1)
                md_pos(k) = new_rnd(Nt,dist_md);
                seq(k,:) = md_pos(k)-(l_md/2-1):md_pos(k)+(l_md/2);
            end
        else
            while any(intersect(seq(1:k-1, :), seq(k,:))) || any(seq(k,:) > Nt) || any(seq(k,:) < 1)
                md_pos(k) = new_rnd(Nt,dist_md);
                seq(k,:) = md_pos(k)-(l_md/2-1):md_pos(k)+(l_md/2);
            end
        end
    end

    for k = 1:n_md
        missing_data(k,:) = signal(seq(k,:));
        missing_data_pos(k,:) = seq(k,:);
    end

end

    function md_pos = new_rnd(Nt,dist_md)
        if strcmp('uniform',dist_md)
            md_pos = randi(Nt,1,1); % uniform distribution
        elseif strcmp('binomial',dist_md)
            md_pos = binornd(Nt,0.5,1,1); % binomial distribution
        else
            error('Unknown distribution')
        end
    end

end
