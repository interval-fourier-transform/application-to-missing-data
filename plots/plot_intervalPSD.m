function p = plot_intervalPSD(w, intervalspectrum)
% Function for visualising the upper and lower bound of the amplitude of
% the interval DFT
%
% INPUT:
%       - w:                Frequency vector
%       - intervalspectrum: Upper and lower bound of the spectrum
%
% OUTPUT:
%       - p:                Plot object with specifications
%
% Marco Behrendt
% Institute for Risk and Reliability, Leibniz Universität Hannover
% behrendt@irz.uni-hannover.de
% https://github.com/marcobehrendt
%
% Date: 17/03/2022

if size(w,1) == 1
    x = [w fliplr(w)];
else
    x = [w; flipud(w)];
end
y = [intervalspectrum(1,:) fliplr(intervalspectrum(end, :))];


face_col_rgb = [0.9020 0.4941 0.3882];
grid_col_rgb = [0.6902 0.6902 0.6902];    

p = fill(x, y, face_col_rgb); 
p.EdgeColor = face_col_rgb; 

grid on;
set(gca, 'GridColor', grid_col_rgb)
set(gca, 'GridAlpha', 1)
set(gca, 'Layer', 'top')

end

