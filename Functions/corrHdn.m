function [ cMax , tMax ] = corrHdn( heading1 , heading2 , begTime , endTime )
% CORRHDN Mean of the dot product between two heading time series. 
% 
%   corrHdn takes two time series of heading, computes the dot product at 
%   every time step, and returns the mean value. This yields a measure, C, 
%   that is analogous to the correlation coefficient, r. To estimate time 
%   delays, C is computed by 'shifting' one heading time series forward or
%   backward in time relative to the other. 
%   
%   Inputs:
%       heading1 = Nx1 vector of heading (in rad); N is # of time steps
%       heading2 = Nx1 vector
%   Outputs:
%       cMax = maximum value of C
%       tMax = shift that produces cMax


% Read heading vectors and decompose them into x,y-components.
x1 = cos(heading1); y1 = sin(heading1);
x2 = cos(heading2); y2 = sin(heading2);
v1 = horzcat(x1,y1);
v2 = horzcat(x2,y2);

% Compute the dot product between v1 and v2 at different 'shifts' to
% determine the effective time delay.
count = 1;
maxShift = 240;
Cij = NaN(2*maxShift+1,2);
for shift = -maxShift:maxShift
    
    % Shift time series; add NaN's to replace data at beginning/end.
    v2s = circshift(v2,shift);
    if shift > 0
        v2s(1:shift,:) = NaN;
    elseif shift < 0
        v2s(end-shift:end,:) = NaN;
    end
    
    % Compute dot product between v1 and shifted v2, w/in specified times.
    v1t = v1(begTime:endTime,:);
    v2t = v2s(begTime:endTime,:);
    dotProductDelay = dot(v1t,v2t,2);
    
    % Store the dot product (effectively, correlation) between v1 and v2,
    % for each shift.
    Cij(count,1) = shift;
    Cij(count,2) = nanmean(dotProductDelay);
    count = count+1;
    
end

% Output maximum dot product (cMax) and the shift that produced it (tMax).
% If cMax is NaN, tMax should also be NaN. 
cMax = max(Cij(:,2));
tMax = Cij(find(Cij(:,2)==cMax),1);
if isnan(cMax) == 1;
    tMax = NaN;
end

end

