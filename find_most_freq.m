function[ans] = find_most_freq(vec)
% This function identifies the most frequently occuring number in vec. If
% more than one number have the highest frequncy, both will be reported.
% ans is a vector of the numbers that have the highest frequency in vec.

nbins = unique(vec);
[a, xout] = hist(vec, nbins);

max_count = max(a);
maxi = find(a == max_count);

ans = xout(maxi);

