function [distance_mat] = pairwise_distances(coordinates)
%pairwise_distances accepts a two-column matrix of x-y co-ordinate values
%(coordinates)and calculates all pair-wise, Euclidean distances between the 
%points specified in coordinates. The output, distance_mat is symmetric, 
%square matrix of pair-wise distances. 

n_points = size(coordinates, 1);
distance_mat = zeros(n_points, n_points) + nan;

%error conditions

%loops for calculating pair-wise distances
for p1n = 1:n_points
    p1 = coordinates(p1n, :);
    for p2n = 1:n_points
        p2 = coordinates(p2n, :);
        dist = sqrt( (p1(1, 1) - p2(1, 1) ).^2 + (p1(1, 2) - p2(1, 2) ).^2   );
        distance_mat(p1n, p2n) = dist;
        clear dist
    end
end