function [AX] = coordAx(point1,point2,point3)
% COORDAX Given 3 points, generates a set of orthogonal unit vectors
% [coordinate axes] centered about the 1st point. The x-axis lies along the
% line connecting point 1 and point 2. The z-axis is orthogonal to the 
% vectors from point 1 to point 2 (v12) and point 1 to point 3 (v13). The 
% y-axis is orthogonal to both the x- and z-axis. 
%   INPUT: 3 points, written in real-world coordinates  
%   OUTPUT: a matrix of 3 vectors, written in real-world coordinates,  
%   which define a marker-centered coordinate system
%
%   EXAMPLE:
%   For the trivial case where the two coordinate systems are the same
%   (i.e. point 1 is located at the origin of the real-world, point 2 is
%   one unit away along the x-axis, and point 3 is one unit away along the
%   y-axis), input would be:
%       point1 = [0 0 0]                                                                          
%       point2 = [1 0 0]
%       point3 = [0 1 0]
%   And output would be
%       AX = [1 0 0; 0 1 0; 0 0 1]
%   Each axis can be called by row:
%       x_axis = AX(1,:) = [1 0 0]
%       y_axis = AX(2,:) = [0 1 0]
%       z_axis = AX(3,:) = [0 0 1]


    % generate vectors connecting x1 to x2 and x3
    v12 = [point2(1)-point1(1) point2(2)-point1(2) point2(3)-point1(3)];    % vector connecting x1 and x2
    v13 = [point3(1)-point1(1) point3(2)-point1(2) point3(3)-point1(3)];    % vector connecting x1 and x3

    % generate 3 orthogonal vectors
    Xnn = v12;                                                              % non-normalized vector defining x-axis in marker-centered coordinate system (written is in real-world coordinates)
    Znn = cross(v12,v13);                                                   % non-normalized vector defining z-axis in marker-centered frame
    Ynn = cross(Znn,Xnn);                                                   % non-normalized vector defining y-axis in marker-centered frame

    % normalize these vectors, to create orthogonal axes of unit length
    X = Xnn/norm(Xnn);                                                      % unit vector along x-axis in marker-centered coordinate system (written in real-world coordinates)
    Y = Ynn/norm(Ynn);                                                      % unit vector along y-axis in marker-centered frame
    Z = Znn/norm(Znn);                                                      % unit vector along z-axis in marker-centered frame    
   
    % combine these vectors into a 3x3 matrix, for output
    AX(:,1) = X;                                                            % first column: vector equation of x-axis (written in real-world coordinates)
    AX(:,2) = Y;                                                            % second column: vector equation of y-axis
    AX(:,3) = Z;                                                            % third column: vector equation of z-axis

    
end

