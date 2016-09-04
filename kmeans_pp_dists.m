function [idx C dists_g dists] = kmeans_pp_dists(X, k)
%This is a function that uses the kmeans ++ algorithm (see function file
%kmeans_pp.m for details) to cluster the m datapoints described in mxn matrix X in n-D 
%space, into k clusters. The outputs idx and C are from the kmeans ++ 
%function. dists_g is a vector of length = no. points with the Euclidean 
%distance of each point from its own cluster-centre in %n-D space. dists is
%a m x k matrix of distances of each point from all cluster centres. 

[idx C] = kmeans_pp(X, k);

no_points = size(X, 1);
no_dims = size(X, 2);
dists_g = zeros(no_points, 1) + nan;
dists = zeros(no_points, k) + nan;

%calculating distances for each point 
for point_no = 1:no_points
    c_num = idx(point_no);              %this point's cluster number
    for c_no = 1:k
        p = X(point_no, :);             %vector describing current point's position
        c = C(:, c_no)';                %vector describing current centre's position
        d = pdist([p; c]);             %distance to current centre
        dists(point_no, c_no) = d;      %saving distance
        
        if c_no == c_num
            dists_g(point_no, 1) = d;   %saving distance to in-cluster distance vec
        else
        end

    end
    
end