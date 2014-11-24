function [x4gb] = lc2gb(x1gb,x2gb,x3gb,x4lc)
% LC2GB Transforms a point x4 in local (marker-centered) coordinates to
% global (real-world) coordinates.
%   INPUT: 3 points (x1,x2,x3), written in local coordinates, and 
%   x4, written in global coordinates  
%   OUTPUT: x4, written in global coordinates
%
%   Use x1,x2,x3 to define the local coordinate system. Then compute a 
%   transformation matrix to convert x4 from local to global coordinates. 
%   The transformation matrix is composed of two matrices: one for rotation
%   and one for translation. 
%
%   Because the two sets of vectors representing the coordinate axes are 
%   orthonormal bases (ONB), computing the rotation matrix turns out to be 
%   relatively easy. As this page explains:
%   http://math.stackexchange.com/questions/23197/finding-a-rotation-transf
%   ormation-from-two-coordinate-frames-in-3-space
%           Since A is a transform from the common "neutral" coordinate 
%       system to the A coordinate system, multiplying with the inverse of 
%       A will be a transform back to the common coordinate system. In this 
%       case, A is an ONB, so "inverse" and "transpose" are synonyms, which 
%       makes your life a lot easier. Similarly, B transforms the neutral
%       coordinate system to the B coordinate system.
%           So, to transform any point from A to B, you multiply it with
%       transpose(A) and then with B.
%
%   The translation matrix is even more straightforward -- it is simply
%   vector addition.

% check if input points are row vectors (1x3) or column vectors (3x1); 
% if column vectors, convert to row vectors
transp = 0;
A = size(x1gb); 
if A(1) == 3
    x1gb = x1gb'; 
    x2gb = x2gb';
    x3gb = x3gb';
end
A = size(x4lc); 
if A(1) == 3
    x4lc = x4lc'; 
    transp = 1; 
end

% generate coordinate axes based on the 3 marker points
AXlc = coordAx(x1gb,x2gb,x3gb);                                             % marker-centered axes                                                
AXgb = [  1  0  0                                                           % real-world axes
          0  1  0
          0  0  1 ];

% rotation matrix
rotMat = transpose(AXlc)*AXgb;

% translation matrix
transMat = x1gb;

% rotate and translate x4 from mc to rw coordinates
x4gb = x4lc*rotMat+transMat;
if transp == 1
    x4gb = x4gb';
end

end

