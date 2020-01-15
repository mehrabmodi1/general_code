function [clust_ids, sorted_ave_traces] = cluster_response_traces(ave_traces, max_n_clusts)
%This function uses matlab's agglomerative clustering algorithm to group
%similar response traces. It then sorts the response traces by cluster
%identity.

%Running clustering code
Z = linkage(ave_traces', 'centroid');
clust_ids = cluster(Z, 'maxclust', max_n_clusts);          %grouping cells into a maximum of 5 clusters
clust_ids_old = clust_ids;
%giving small clusts the highest clust-id numbers
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


%re-ordering responses acc to cluster identity
sorted_resps = [];
for clust_n = 1:n_clusts
    curr_cells = find(clust_ids == clust_n);
    curr_resps = ave_traces(:, curr_cells);       %response vectors for currently clustered cells
    ave_resp = mean(curr_resps, 2, 'omitnan');

    corrs = corrcoef([ave_resp, curr_resps]);
    corrs = corrs(1, 2:end);
    curr_resps = [corrs; curr_resps];
    curr_resps = sortrows(curr_resps', -1)';
    sorted_resps = [sorted_resps, curr_resps(2:end, :)];

end
keyboard