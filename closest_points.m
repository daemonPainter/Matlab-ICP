function [ Y ] = closest_points( P, X )
% Find closest points between two sets of points.
%   P, X must be 4-by-n in homogeneous coordinates or 3-by-n in cartesian
%   coordinates.
%   [Y] = closest_points(P,X) will reorder points of P so they match the
%   same order of the relative closest neighbour in X and will save the
%   re-ordered P set in Y.
%   Works only if X and P have the same number of points. See issue #3 on
%   repository if you'd like to contribute.

P = P(1:3,:);
X = X(1:3,:);

n_points = size(P,2);
Y = zeros(3,n_points);

% mean value to remove errors and disturbances
P_mean = mean(P,2);
X_mean = mean(X,2);

P = P - repmat(P_mean, 1, n_points);
X = X - repmat(X_mean, 1, n_points);

for i = 1:n_points
    % create a matrix the same size of P with only one point of X,
    diff_matrix = (repmat(X(:,i), 1, n_points) - P);
    % find distances between all these points,
    distance = sqrt(sum(diff_matrix.^2, 1));
    % we are interested in the closest one,
    [~, index] = min(distance);
    % select that points (add mean value)
    Y(:, i) = P(:,index) + P_mean;
end
