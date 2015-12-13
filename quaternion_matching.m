function [T] = quaternion_matching(P,X)
%Quaternion Matching computes homogeneous transformation
%   [T] = quaternion_matching(P,X)
%   P is a measured data point set to be aligned with a model point set X.
%   P and X must have homogeneous coordinates and be 4-by-n
%   (e.g. [x;y;z;1]), where n is the number of points. T is the transform
%   in homogeneous coordinates that will register points in P so they are
%   expressed in the same reference frame as X. If P and X are the same set
%   of points expressed in two different reference frames, then X = T*P.

P = P(1:3,:);   % back to cartesian coordinates
X = X(1:3,:);

N_P = size(P,2);    % number of points in P
N_X = size(X,2);    % number of points in X, currently unused

ni_P = mean(P,2);    % Center of mass
ni_X = mean(X,2);

% Cross-covariance matrix
Sigma_px = zeros(3);
for i = 1 : N_P
    Sigma_px = Sigma_px + ((P(:,i) - ni_P) * (X(:,i) - ni_X)');
end
Sigma_px = 1 / N_P * Sigma_px;

% Compute parameters for the quaternion extraction
    % Anti-symmetric matrix
    A = Sigma_px - Sigma_px';
    % cyclic components of the anti-symmetric matrix
    D = [A(2,3) A(3,1) A(1,2)]';
    % 4-by-4 symmetric matrix Q(Sigma_px)
    Q = [trace(Sigma_px) D';...
         D               (Sigma_px + Sigma_px' - trace(Sigma_px) * eye(3))];

% the unit eigenvector q_R corresponding to the maximum eigenvalue of the
% matrix Q(Sigma_px) is selected as the optimal rotation
    % Eigenval/Eigenvect extraction
    [E_vec,E_val] = eig(Q);

    % Find the maximum eigenval and its position
    [~,index] = max(diag(E_val));

    % Extract q_R
    q_R = E_vec(:,index);
    
% Calculates the quaternion <- this is a class
quat = Quaternion(q_R);

% Extract the homogeneous transform <- class method
T = quat.T; % <- has no trans component

% The optimal translation vector is given by
T(1:3,4) = (ni_X - (quat.R) * ni_P)';   % <- using class method
