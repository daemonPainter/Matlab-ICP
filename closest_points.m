function [ cp ] = closest_points( source, model )
% Find closest points between two sets of points.
%   Input must be 4xN. Will reorder points of source so they match the same
%   index of the closest neighbour in model.

source = source(1:3,:);
model = model(1:3,:);

n_points = size(source,2);

source_mean = mean(source,2);

source = source - repmat(source_mean, 1, n_points);
model = model - repmat(mean(model, 2), 1, n_points);

for i = 1:n_points
    diff_matrix = (repmat(model(:,i), 1, n_points) - source);
    distance = sqrt(sum(diff_matrix.^2, 1));
    [~, index] = min(distance);
    cp(:, i) = source(:,index) + source_mean;

end

% calcolo media dei punti dei due set -> avvicino lungo le dimensioni i
% punti corrispondenti mean(data,2) repmat su tutti gli altri 
% sistema può funzionare anche con numero diverso tra modello e dati.
% media dei punti sulle tre coordinate  e poi aggiungo quelle della source
% nell'output
