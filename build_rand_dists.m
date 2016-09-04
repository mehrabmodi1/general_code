function [dist_mat, count_mat] = build_rand_dists(frac_vec, n_cells, n_dists)
 
% frac_vec = [0.2, 0.1, 0.4];
% n_cells = 100;
% n_dists = 1000;

    n_vars = length(frac_vec);                      %the variables for which independent counts are known, such as odors 
    dist_mat = zeros(n_cells, n_vars, n_dists);     %the matrix that will contain randomly drawn cell populations, positive for each variable
    rand_mat = rand(n_cells, n_vars, n_dists);      %the matrix of random numbers used to generate dist_mat

    
    %loop to randomly pick responder cells for each variable
    for var_n = 1:n_vars
        curr_frac = frac_vec(var_n);
        a = find(rand_mat(:, var_n, :) < curr_frac);                    %find a random set of cells that on average, are positive for the current var at curr_frac proportion
        [cell_n, mat_n] = ind2sub(size(rand_mat(:, var_n, :)), a);      
        
        for dist_mat_n = 1:n_dists
            curr_cellsi = find(mat_n == dist_mat_n);
            dist_mat(cell_n(curr_cellsi), var_n, dist_mat_n) = 1;
        end
        
    end
    
    var_pairs = combinator(n_vars, 2, 'c');
    %counting cell-fractions for combinations of vars
    count_mat = zeros(size(var_pairs, 1), n_dists) + nan;              %matrix to store counts of cells responding to each pair of stims
    for pair_n = 1:size(var_pairs, 1)
        curr_pair = var_pairs(pair_n, :);
        sum_vecs = dist_mat(:, curr_pair(1), :) + dist_mat(:, curr_pair(2), :);
        
        both_positives = find(sum_vecs == 2);
        [cell_n, mat_n] = ind2sub(size(sum_vecs), both_positives);
        for dist_mat_n = 1:n_dists
            both_positivesi = find(mat_n == dist_mat_n);
            n_both_positives = size(both_positivesi, 1);
            count_mat(pair_n, dist_mat_n) = n_both_positives;
        end
        
    end

end