function [corr_val] = mat_corrcoef(mat1, mat2)
%this fucntion reshapes 2D matrices into single vectors and then
%calculates the corrcoef of the two vectors

vec1 = reshape(mat1, 1, []);
vec2 = reshape(mat2, 1, []);

corr_val = corrcoef(vec1, vec2);
corr_val = corr_val(1, 2);

end