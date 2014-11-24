function [ speed , heading ] = xyToSpeedHeading(x,y)
% XYTOSPEEDHEADING Produce speed and heading time series from time series 
% of position (x,y). 

speed = sqrt((diff(x).^2+diff(y).^2))/(1/60);
heading = atan2(diff(x),diff(y));

speed = vertcat(NaN,speed);  
heading = vertcat(NaN,heading);

end