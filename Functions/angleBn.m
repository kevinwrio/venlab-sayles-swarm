function [ angle_bn_angles ] = angle_bn(a,b)
% ANGLE_BN_ANGLES computes the angle between two angles.
%   Equation borrowed from: 
%   http://www.mathworks.com/matlabcentral/newsreader/view_thread/151925

angle_bn_angles = atan2(sin(a-b),cos(a-b));

end

