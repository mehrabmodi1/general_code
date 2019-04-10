function [corr] = mat_corrcoef(mat1, mat2)
% this function reshapes two 2d matrices into long vectors and calculates
% the correlation coefficients of the two vectors


vec1 = reshape(mat1, [], 1);
vec2 = reshape(mat2, [], 1);

c = corrcoef(vec1, vec2);
corr = c(1, 2);


end
