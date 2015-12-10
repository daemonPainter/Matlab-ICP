function [T] = quaternion_matching(Pa,Pb)
%Quaternion Matching computes homogeneous transformation T so that Pb=T*Pa
% Pa, Pb must be homogeneous coordinates with dimension 4 x n, where n is
% the number of points

Pa = Pa(1:3,:);   % Remove 1s
Pb = Pb(1:3,:);

dims = size(Pa);
n_points = dims(2);  % Number of points

Ba = mean(Pa,2);    % Center of mass
Bb = mean(Pb,2);

% Cross correlation function
cross_matrix = zeros(dims(1));
for i = 1 : n_points
    cross_matrix = cross_matrix + ((Pa(:,i) - Ba) * (Pb(:,i) - Bb)');
end
cross_matrix = 1 / n_points * cross_matrix;

% Calculate the parameters for the quaternion extraction
A = cross_matrix - cross_matrix';
D = [A(2,3) A(3,1) A(1,2)]';
Q = [trace(cross_matrix), D'; D, (cross_matrix + cross_matrix' - trace(cross_matrix) * eye(3))];

% Eigenval/Eigenvect extraction
[E_vec,E_val] = eig(Q);

% Find the maximum eigenval and its position
[~,max_p] = max(diag(E_val));

% Calculates the quaternion
quat = Quaternion(E_vec(:,max_p));

% Extract the rotation matrix
T = quat.T;

% Calculates the translation component
T(1:3,4) = (Bb - (quat.R) * Ba)';
