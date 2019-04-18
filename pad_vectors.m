function [mat_out] = pad_vectors(vectors, nanpad);
%This function takes a cell-array as input, with vectors of different
%lengths in different cells (one in each cell). It then identifies the
%longest vector, pads all the others with NaNs (if nan == 1) or 0's if nan
%== 0 and concatenates the vectors into a single matrix, with n columns
%equal to n vectors. NOTE: vector cells must be arranged in a single row in
%'vectors'. All vectors must be columns.
%syntax: mat_out = pad_vectors(vectors, nan);


n_vecs = max(size(vectors));     %number of vectors to pad;

%finding out length of longest vector
lengths = [];
for vec_n = 1:n_vecs
    curr_vec = vectors{1, vec_n};
    lengths = [lengths; length(curr_vec)];
    
end
l_adjust = max(lengths);        %length of longest vector


%saving vectors to mat_out
mat_out = zeros(l_adjust, n_vecs);
if nanpad == 1
    mat_out = mat_out + nan;
else
end

for vec_n = 1:n_vecs
    curr_vec = vectors{1, vec_n};
    curr_l = length(curr_vec);
    mat_out(1:curr_l, vec_n) = curr_vec;
    
end

