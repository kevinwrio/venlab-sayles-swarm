function [vOutput] = bFiltExtrap (hZ,vInput)
% BFILTEXTRAP   Butterworth low-pass filter with linear extrapolation. 
%   butterFilt(hZ, varargin) Filters and returns input vector with a 
%   Butterworth low-pass filter using standard VENLab parameters and 
%   sample rate hZ. Extrapolates last 0.5 sec to avoid transient at end. 

% Created by Henry Harrison, 2010
% Modified by Kevin Rio, 2013
% Brown University

% Compute parameters
% cutoff = (0.6/0.802);                                                       % standard VENLab parameters
cutoff = 1;                                                                 % parameters from Rio, Rhea, & Warren (2014)   
[B,A] = butter(4,cutoff/(hZ/2));                                      

% Determine how many time-steps to use for extrapolation (default = 30)
if length(vInput) > 30
    pad = 30;
else
    pad = length(vInput);
end

% Extrapolation
vLine = linspace(1,pad,pad)';
vEnd = vInput(end-(pad-1):end,1);                                           % end of input vector // extrapolate from this
eqn = polyfit(vLine,vEnd,1); m = eqn(1); b = eqn(2);                        % equation of best fit line (y = mx+b) to vEnd
ext = 1:pad*4;
vExtension = (m*ext+b+(vInput(end)-vInput(end-(pad-1))))';                  % data points created by extrapolation
vInExtended = vertcat(vInput,vExtension);                                   % input vector, extended with these data points

% Apply filter and return input to original size
vOutExtended = filtfilt(B,A,vInExtended);                                   % filtered time series will be longer than input ...
vOutput = vOutExtended(1:end-pad*4);                                        % ... so truncate back to original length

