function T = ECE569_MatrixExp6(se3mat)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes a se(3) representation of exponential coordinates.
% Returns a T matrix in SE(3) that is achieved by traveling along/about the 
% screw axis S for a distance theta from an initial configuration T = I.
% Example Input:
% 
% clear; clc;
% se3mat = [ 0,      0,       0,      0;
%          0,      0, -1.5708, 2.3562;
%          0, 1.5708,       0, 2.3562;
%          0,      0,       0,      0]
% T = MatrixExp6(se3mat)
% 
% Output:
% T =
%    1.0000         0         0         0
%         0    0.0000   -1.0000   -0.0000
%         0    1.0000    0.0000    3.0000
%         0         0         0    1.0000 

omgtheta = ECE569_so3ToVec(se3mat(1:3, 1:3));
if ECE569_NearZero(norm(omgtheta))
    T = [eye(3), se3mat(1:3, 4); 0, 0, 0, 1];
else
    [~, theta] = ECE569_AxisAng3(omgtheta);
    omgmat = se3mat(1:3, 1:3) / theta;
    R = ECE569_MatrixExp3(se3mat(1:3, 1:3));
    v = se3mat(1:3, 4);
    G = (eye(3)*theta + (1-cos(theta))*omgmat + (theta-sin(theta))*omgmat*omgmat) / theta;
    T = [R, G*v; 0, 0, 0, 1];
end
end