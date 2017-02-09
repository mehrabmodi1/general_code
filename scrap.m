n_cells = 10000;
n_odors = 100;
sparseness = .1;

sig_resp_mat = (rand(n_cells, n_odors));

del = find(sig_resp_mat >= (1 - sparseness));
sig_resp_mat(del) = 1;
del = find(sig_resp_mat < (1 - sparseness) );
sig_resp_mat(del) = 0;

figure(1)
a = histogram(sum(sig_resp_mat, 2))

% figure(2)
% b = histogram(sum(sig_resp_mat, 1))