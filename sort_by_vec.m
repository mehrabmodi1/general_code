function [mat_out] = sort_by_vec(mat, vec, asc)
%Syntax: [mat_out] = sort_by_vec(mat_in, vec, asc)
%This function sorts mat_in by the entries in vec. Asc can be 1 or 0 and 
%specifies whether the sorting is ascending(1) or descending (0). 

siz_vec = size(vec);

if siz_vec(2) > siz_vec(1)
    vec = vec';
    siz_vec = size(vec);
else
end

if size(mat, 2) == size(vec, 1)
    flipped_mat = 1;
    mat = mat';
else
end

mat = [vec, mat];

if asc == 0
    mat = sortrows(mat, 'descend');
else
    mat = sortrows(mat, 'ascend');
end

if flipped_mat == 1
    mat = mat';
else
end

mat_out = mat(:, 2:size(mat, 2));
    