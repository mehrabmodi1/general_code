function [grouped_vals_vec2] = bin_vec2_by_vec1(refvec, vec2, bins)
%Syntax: function [vec2_binned] = bin_vec2_by_vec1(refvec, vec2, bins)
%This function uses the values in vec1 as a reference, and bins them into
%pre-specified bins. It then identifies the the values in vec2 at the same 
%indices, and puts them in a corressponding set of bins in vec2 space. The
%values in column n in grouped_vals_vec2 consist all the values in vec2
%that have the same indices as values in refvec that fall in bin n.

%testing lines
% refvec = 1:100;
% vec2 = [ones(1, 20), zeros(1, 50), ones(1, 30)];
% bins = 0:10:100;

%standarising vec dimensions
if size(refvec, 1) == 1
    refvec = refvec';
else
end
if size(vec2, 1) == 1
    vec2 = vec2';
else
end

%grouping vec2 values according to indices of refvec that belong in each bin
grouped_vals_vec2 = [];
for bin_n = 2:length(bins)
    curr_i = find(refvec > bins(bin_n - 1) & refvec < bins(bin_n));    %indices in ref_vec that fall in current bin
    if isempty(curr_i) == 0
        grouped_vals_vec2 = pad_n_concatenate(grouped_vals_vec2, vec2(curr_i), 2, nan);
    else
        grouped_vals_vec2 = pad_n_concatenate(grouped_vals_vec2, nan, 2, nan);
    end
    
end


