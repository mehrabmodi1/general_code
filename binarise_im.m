function [bin_im] = binarise_im(im, thresh)
%this function is redundant to im2bw, except that it uses an absolute value
%threshold to binarise images, rather than a threshold scaled between 0 and
%1.

bin_im = zeros(size(im));
a = find(im > thresh);
bin_im(a) = 1;


end

