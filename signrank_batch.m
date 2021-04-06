function [p_vals] = signrank_batch(data_mat)
%syntax: [p_vals] = signrank_batch(data_mat)
%This function calls signrank in a loop to do a Wilcoxon signed rank test
%for each column vector in data_mat to check if it's median is different
%from 0
p_vals = zeros(1, size(data_mat, 2));
for vec_n = 1:size(data_mat, 2)
    if sum(isnan(data_mat(:, vec_n))) == size(data_mat, 1)
        continue
    else
    end
    p_vals(1, vec_n) = signrank(data_mat(:, vec_n));    
end