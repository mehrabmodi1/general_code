function [concatenated_mat] = concatenate_padded(mat1, mat2, concatenate_dim, pad_with)
%syntax: [concatenated_mat] = concatenate_padded(mat1, mat2, concatenate_dim, pad_with)
%This function concatenates two matrices of unequal sizes along the
%dimension concatenate_dim, padding with the number specified in pad_with (typically nan or 0).

% clear all
% close all
% 
% mat1 = [1, 1, 1; 1, 2, 1];
% mat2 = zeros(3, 3, 3) + 1;
% concatenate_dim = 2;
% pad_with = nan;

size1 = size(mat1);
size2 = size(mat2);

n_common_dims = min([length(size1), length(size2)]);
for dim_n = 1:n_common_dims
    longer_dim = max([size1(dim_n), size2(dim_n)]);
    dim_lengths(1, dim_n) = longer_dim;
    
end

n_extra_dims = abs(length(size1) - length(size2));
larger_mat = sign(length(size1) - length(size2));

if abs(larger_mat) > 0
    if larger_mat == 1
        size_x = size1;
    elseif larger_mat == -1
        size_x = size2;
    else
    end
    
    for dim_n = (n_common_dims + 1):(n_common_dims + n_extra_dims)
        dim_lengths = [dim_lengths, size_x(dim_n)];
        
    end
else
end

%creating an empty matrix the size of the final concatenated matrix and
%replacing 0's with pad_with
dim_lengths(1, concatenate_dim) = size1(concatenate_dim) + size2(concatenate_dim);
concatenated_mat = zeros(dim_lengths) + pad_with;

%putting mat1 into padded matrix
%assmebling execution string
exec_string = [];
for dim_n = 1:length(size1)
    if dim_n > 1
        exec_string = [exec_string ','];
    else
    end
    exec_string = [exec_string '1:(' int2str(size1(dim_n)) ')'];
end
eval(['concatenated_mat(' exec_string ') = mat1;' ]);


%putting mat2 into padded matrix, after mat1 along concatenate_dim
%assmebling execution string
exec_string = [];
for dim_n = 1:length(size2)
    if dim_n > 1
        exec_string = [exec_string ','];
    else
    end
    if dim_n == concatenate_dim
        exec_string = [exec_string int2str(size1(concatenate_dim) + 1) ':' int2str(size(concatenated_mat, dim_n))];
    elseif dim_n ~= concatenate_dim
        exec_string = [exec_string '1:' int2str(size2(dim_n))];
    else
    end
end
eval(['concatenated_mat(' exec_string ') = mat2;' ]);




%end