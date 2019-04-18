function [com_vec] = compute_COM_vec(response_mat, t_dim)
%syntax: [com_vec] = compute_COM_vec(response_mat, t_dim)
%This function loops through all the cells in response_mat and returns a
%vector of the centers of mass of the response vectors in time. t_dim is
%the dimension number of time in response_mat.

if t_dim == 2
    response_mat = response_mat';
else
end

n_cells = size(response_mat, 2);

com_vec = zeros(1, n_cells) + nan;
for cell_n = 1:n_cells
    resp_vec = response_mat(:, cell_n);
    min_val = min(resp_vec);
    if sign(min_val) == -1
        resp_vec = resp_vec + abs(min_val).*1.001 ;
        resp_vec(isnan(resp_vec)) = [];
    else
    end
    try
        com = centerOfMass(resp_vec);
    catch
        keyboard
    end
    com_vec(1, cell_n) = com(1, 1);
    
end