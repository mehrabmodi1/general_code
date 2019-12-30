function [wtd_mean] = weighted_mean(x, wts)
%This function computes the mean of each column in x, where each row is
%weighted by the weights specified in the vector wts. length(wts) should
%equal n rows in x.

n_cols = size(x, 2);
x_wt = x.*repmat(wts, 1, n_cols);

wtd_mean = sum(x_wt, 1);
wtd_mean = wtd_mean./sum(wts);