function [b1,phase,N] = rdpls(fname)
% RDPLS - reads a pulse file in VNMR format, converts the values to 
%	the form expected by PJB's Matlab functions, described below.
%
% Expected input format:
% 	There should be three columns, delimted by tabs (or space):
%		phase (degrees), B1 (unitless), and freq(ignored)
% 	This is the same format used by Varian's pulsetool (perhaps VNMR too?)
%	Pulsetool requires that B1 have a max of 1024. This function, however, 
%		renormalizes B1 for its output. 
% 
% Output:
% 	b1 - pulse amplitude, normalized to 1. 
% 	phase - phase, converted to radians
% 	N - number of steps in pulse

% PJB 06.10.00
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Raw ascii read
pulse = load(strcat(fname));

% The array size is the number of lines of input
N=size(pulse,1);

% Renormalize B1 to 1. Unitless.
b1max = max(pulse(:,2));
b1 = pulse(:,2)/(1.0*b1max);

% Phase is input in degrees, output in radians.
phase = pulse(:,1)*pi/180;