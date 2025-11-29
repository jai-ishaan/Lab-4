function expmat = ECE569_MatrixLog6(T)
% *** CHAPTER 3: RIGID-BODY MOTIONS ***
% Takes a transformation matrix T in SE(3).
% Returns the corresponding se(3) representation of exponential 
% coordinates.
% Example Input:
% 
% clear; clc;
% T = [[1, 0, 0, 0]; [0, 0, -1, 0]; [0, 1, 0, 3]; [0, 0, 0, 1]];
% expmat = MatrixLog6(T)
% 
% Output:
% expc6 =
%         0         0         0         0
%         0         0   -1.5708    2.3562
%         0    1.5708         0    2.3562
%         0         0         0         0

[R, p] = ECE569_TransToRp(T);
omgmat = ECE569_MatrixLog3(R);
if isequal(omgmat, zeros(3))
    expmat = [zeros(3), p; 0, 0, 0, 0];
else
    theta = acos((trace(R)-1)/2);
    omgmat_unnorm = omgmat / theta;
    a = theta / 2;
    b = (1 - cos(theta)) / theta;
    c = (theta - sin(theta)) / theta;
    G_inv = eye(3) * a + omgmat_unnorm * b + omgmat_unnorm * omgmat_unnorm * c;
    G_inv = G_inv * theta;
    v = G_inv * p;
    expmat = [omgmat, v; 0, 0, 0, 0];
end
end