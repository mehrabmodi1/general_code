function [overlap_areas] = dist_overlap_calc(bins, dist1, dist2)
%This function calculates the area of overlap between two arbitrary
%distributions, calculated for the same vector of bins. The output is the
%overlapping area at each bin value. Total overlapping area is the sum of
%this vector.
%Syntax:[overlap_areas] = dist_overlap_calc(bins, dist1, dist2)

if length(bins) == length(dist1) && (length(bins) + 1) == (length(dist2) + 1)
else
    error('number of bin edges should be length of each dist plus 1.')
end

n_bins = length(bins - 1);
overlap_areas = zeros(n_bins, 1);
for bin_n = 1:n_bins
    overlap_areas(bin_n) = min([dist1(bin_n), dist2(bin_n)]);
end