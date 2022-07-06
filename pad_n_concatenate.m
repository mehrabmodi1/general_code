function [conc_mat] = pad_n_concatenate(big_mat, small_mat, dim_conc, pad_with)
%syntax:  [conc_mat] = pad_n_concatenate(big_mat, small_mat, dim_conc, pad_with)
%This function pads big_mat or small_mat to match them in size in all
%dimensions except dim and then concatenates small_mat to the end of big_mat
%along the dimension dim_conc. The pads consist of extra rows/columns filled
%with pad_with (typically 0 or NaN).
%Mehrab Modi 20160617
 
% big_mat = [];
% small_mat = [5, 5; 5, 5];
% dim_conc = 3;
% pad_with = nan;

size_big = size(big_mat);
size_small = size(small_mat);

n_dims = max([size(size_big, 2), size(size_small, 2), dim_conc]);       %the number of dimensions conc_mat will eventually have

%checking if either matrix has a smaller number of dimensions and padding
%size vecs appropriately
dif_big = n_dims - size(size_big, 2);
if dif_big > 0
    size_big = [size_big, zeros(1, dif_big)];
else
end
dif_small = n_dims - size(size_small, 2);
if dif_small > 0
    size_small = [size_small, zeros(1, dif_small)];
else
end

max_lengths = max([size_big; size_small]);          %picking out the longest lengths along each dimension between big_mat and small_mat



%Initialising the final conc_mat as an empty matrix
non_zeros = find(max_lengths > 0);      %counting the number of non-zero dims in max_lengths
if length(non_zeros) < dim_conc
    if max(size_big) > 0
        max_lengths(dim_conc) = 2;
    elseif max(size_big) == 0
        max_lengths(dim_conc) = 1;
    else
    end
    conc_mat_siz = max_lengths;
else
    conc_mat_siz = max_lengths;
    if size_small(dim_conc) > 0
        conc_mat_siz(dim_conc) = size_big(dim_conc) + size_small(dim_conc);
    elseif size_small(dim_conc) == 0
        conc_mat_siz(dim_conc) = size_big(dim_conc) + 1;
    else
    end
end

try
    conc_mat = zeros(conc_mat_siz) + pad_with;
    conc_mat_init = conc_mat;
catch
    %if you're getting an error here, make sure you're using nan and not
    %'nan' as pad_with
end

%filling out conc_mat with big_mat and small_mat, with big_mat going first
%along dim_conc and small_mat going second

%building eval string to fill in big_mat
big_mat_string = ['conc_mat('];
small_mat_string = ['conc_mat('];
for dim_n = 1:n_dims
    if dim_n == dim_conc
        if size_big(dim_n) > 0
            big_mat_string = [big_mat_string, '1:',  int2str(size_big(dim_n)), ','];
        elseif size_big(dim_n) == 0
            big_mat_string = [big_mat_string, '1:1,'];
        end
        
        if size_small(dim_n) > 0
            small_mat_string = [small_mat_string, int2str(size_big(dim_n) + 1), ':' , int2str(size(conc_mat, dim_n)), ','];
        elseif size_small(dim_n) == 0 
            if max(size_big) > 0
                if size_big(dim_n) == 0
                    small_mat_string = [small_mat_string, '2:2,'];
                elseif size_big(dim_n) > 0
                    small_mat_string = [small_mat_string, int2str(size_big(dim_n) + 1) ':' int2str(size_big(dim_n) + 1), ','];
                else
                end
            else
                small_mat_string = [small_mat_string, '1:1,'];
            end
        end
        
    elseif dim_n ~= dim_conc
        if size_big(dim_n) > 0
            big_mat_string = [big_mat_string, '1:',  int2str(size_big(dim_n)), ','];
        else
        end
        
        if size_small(dim_n) > 0
            small_mat_string = [small_mat_string, '1:',  int2str(size_small(dim_n)), ','];
        else
        end
        
    else
    end
    
end

big_mat_string = [big_mat_string(1:(end-1)), ') = big_mat;'];
small_mat_string = [small_mat_string(1:(end-1)), ') = small_mat;'];

if max(size_big) > 0
    eval(big_mat_string)
else
end

if max(size_small) > 0
    eval(small_mat_string)
else
end


end