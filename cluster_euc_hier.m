function [data_mat_sorted, data_mat_sorted_n, label_mat_sorted, clust_ids, link_map] = cluster_euc_hier(data_mat, label_mat, n_clusts)
%Syntax:function [data_mat_sorted, data_mat_sorted_n, label_mat_sorted, clust_ids, link_map] = cluster_euc_hier(data_mat, label_mat, n_clusts)
%This function takes a data matrix with nrows observations and ncols
%dimensions and clusters observations by their euclidean distances with the 
%linkage function. n_clusts is the maximum number of clusters expected.
%It then assigns clusters with the cluster function, and returns a sorted 
%corrcoef matrix for display, along with cluster ids and the linkage map.
%labelmap is a matrix of nrows = nrows in data_mat. It is an arbitrary set
%of labels to be co-sorted along with the observations.

%normalising columns of data_mat by computing z-scores
mean_vec = mean(data_mat, 'omitnan');
sd_vec = std(data_mat, [], 'omitnan');
data_mat_orig = data_mat;
data_mat = data_mat - repmat(mean_vec, size(data_mat, 1), 1);
data_mat = data_mat ./ repmat(sd_vec, size(data_mat, 1), 1);

del = isnan(data_mat);
data_mat(del) = 0;

if isempty(label_mat) == 1
    label_mat = zeros(size(data_mat, 1), size(data_mat, 2));
else
end

link_map = linkage(data_mat, 'centroid');
clust_ids = cluster(link_map, 'maxclust', n_clusts);          %grouping cells into a maximum of 5 clusters
clust_ids_old = clust_ids;

%giving small clusts the largest clust-id numbers
n_clusts = max(clust_ids);
n_cells_list = zeros(n_clusts, 2);
for clust_n = 1:n_clusts
    n_cells_list(clust_n, 2) = clust_n;
    n_cells_list(clust_n, 1) = length(find(clust_ids == clust_n));
end
n_cells_list = sortrows(n_cells_list, -1);

for clust_n = 1:n_clusts
    old_clust_n = n_cells_list(clust_n, 2);
    curr_cells = find(clust_ids_old == old_clust_n);
    clust_ids(curr_cells) = clust_n;
end


%re-ordering observations acc to cluster identity
data_mat_sorted = [];
data_mat_sorted_n = [];
label_mat_sorted = [];
clust_ids_sorted = [];
for clust_n = 1:n_clusts
    curr_members = find(clust_ids == clust_n);
    curr_labels = label_mat(curr_members, :);   %label vectors for currently cluster's observations
    curr_obs = data_mat_orig(curr_members, :);       %observation vectors for currently cluster's observations
    curr_obs_n = data_mat(curr_members, :);       %observation vectors for currently cluster's observations
    ave_obs = mean(curr_obs, 1, 'omitnan');    
    
    %sorting within current cluster
    try
        corrs = corrcoef([ave_obs; curr_obs]');
    catch
        keyboard
    end
    corrs = corrs(2:end, 1);
   
    curr_obs = [corrs, curr_obs];
    curr_obs = sortrows(curr_obs, -1);
    data_mat_sorted = [data_mat_sorted; curr_obs(:, 2:end)];
    
    curr_labels = [corrs, curr_labels];
    curr_labels = sortrows(curr_labels, -1);
    label_mat_sorted = [label_mat_sorted; curr_labels(:, 2:end)];
    
    curr_obs_n = [corrs, curr_obs_n];
    curr_obs_n = sortrows(curr_obs_n, -1);
    data_mat_sorted_n = [data_mat_sorted_n; curr_obs_n(:, 2:end)];
        
    clust_ids_sorted = [clust_ids_sorted; (zeros(length(curr_members), 1) + clust_n) ];
end
%clust_ids = clust_ids_sorted;