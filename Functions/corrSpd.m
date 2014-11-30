function [ rMax , tMax ] = corrSpd( speed1 , speed2 , begTime , endTime )
% CORRSPD Cross-correlation between two heading time series. 
% 
%   corrHdn takes two time series of heading, and returns the
%   cross-correlation (coefficient, r, and delay, t). To estimate time
%   delays, the time series are 'shifted' forward or backward in time
%   relative to one another. 
%   
%   Inputs:
%       speed1 = Nx1 vector of heading (in rad); N is # of time steps
%       speed2 = Nx1 vector
%   Outputs:
%       rMax = maximum value of r
%       tMax = shift that produces rMax


% Read and rename speed vectors.
v1 = speed1;
v2 = speed2; 

% Compute the dot product between v1 and v2 at different 'shifts' to
% determine the effective time delay.
count = 1;
maxShift = 240;
Rij = NaN(2*maxShift+1,2);
for shift = -maxShift:maxShift
    
    % Shift time series; add NaN's to replace data at beginning/end.
    v2s = circshift(v2,shift);
    if shift > 0
        v2s(1:shift) = NaN;
    elseif shift < 0
        v2s(end-shift:end) = NaN;
    end
    
    % Compute dot product between v1 and shifted v2, w/in specified times.
    v1t = v1(begTime:endTime);
    v2t = v2s(begTime:endTime);
    
    % Compute correlation between v1 and shifted v2, w/in specified times.
    corrMatrix = corrcoef(v1t,v2t,'rows','complete');
    corrDelay = corrMatrix(2);
    
    % Store the correlation between v1 and v2, for each shift. 
    Rij(count,1) = shift;    
    Rij(count,2) = corrDelay; 
    count = count+1;
    
end

% Output maximum correlation (rMax) and the shift that produced it (tMax).
rMax = max(Rij(:,2));
tMax = Rij(find(Rij(:,2)==rMax),1);
if isnan(rMax) == 1;
    tMax = NaN;
end

end

