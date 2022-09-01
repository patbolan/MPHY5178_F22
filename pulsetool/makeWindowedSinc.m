% Makes a windowed sinc waveform
%
% INPUTS
%   Npts: number of points in waveform
%   Nlobes: number of lobes in sinc. Typically odd (3,5,7) 
%   alpha: window multiplier
%       1 = no window
%       0.5 = Hanning Window
%       0.54 = Hamming Window
%
% OUPUTS
%   b1 = the waveform
function b1 = makeWindowedSinc(Npts, Nlobes, alpha)

% tau goes from -0.5 to 0.5, with max at tau=0
tau = ((0:Npts-1) ./ (Npts-1)) - 0.5;

% Make the sinc function
b1 = sin((Nlobes+1) * pi .* tau) ./ ((Nlobes+1) * pi .* tau);

% at tau==1 the denominator is NaN, b1 should be 1
b1(isnan(b1)) = 1;

% Calculate and apply window function
window = alpha + (1-alpha)*cos(2 * pi .* tau);
b1 = b1.*window;

