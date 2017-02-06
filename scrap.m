n_cells = 1000;
n_odors = 5;
sparseness = .1;

sig_resp_mat = (rand(10000, 5));

del = find(sig_resp_mat >= 0.9);
sig_resp_mat(del) = 1;
del = find(sig_resp_mat < 0.9);
sig_resp_mat(del) = 0;

histogram(sum(sig_resp_mat, 2))