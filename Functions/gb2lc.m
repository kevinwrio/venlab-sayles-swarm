function [x4lc] = gb2lc(x1gb,x2gb,x3gb,x4gb)
% TRANSMAT Transforms a point x4 in global (real-world) coordinates to
% local (marker-centered) coordinates. 
%   INPUT: 3 points (x1,x2,x3), written in global coordinates, and 
%   x4, written in global coordinates  
%   OUTPUT: x4, written in local coordinates
%
%   Use x1,x2,x3 to define the local coordinate system. Then compute a 
%   transformation matrix to convert x4 from global to local coordinates. 
%   See lc2gb for mathematical details. 

% check if input points are row vectors (1x3) or column vectors (3x1); 
% if column vectors, convert to row vectors
transp = 0;
A = size(x1gb); 
if A(1) == 3
    x1gb = x1gb'; 
    x2gb = x2gb';
    x3gb = x3gb';
end
A = size(x4gb); 
if A(1) == 3
    x4gb = x4gb'; 
    transp = 1; 
end

% generate coordinate axes based on the 3 marker points
AXlc = coordAx(x1gb,x2gb,x3gb);                                             % marker-centered axes                                                
AXgb = [  1  0  0                                                           % real-world axes
          0  1  0
          0  0  1 ];

% rotation matrix
rotMat = transpose(AXgb)*AXlc;

% translation matrix
transMat = -x1gb;

% rotate and translate x4 from mc to rw coordinates
x4lc = (x4gb+transMat)*rotMat;
if transp == 1
    x4lc = x4lc';
end

end

