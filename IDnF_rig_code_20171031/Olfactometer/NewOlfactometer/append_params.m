function [params_mat_fused] = append_params(params_mat1, params_mat2)

n_trs1 = size(params_mat1.duration, 1);
n_trs2 = size(params_mat2.duration, 1);

params_mat_fused = params_mat1;

%building list of fields to append
field_lengths = structfun(@numel,params_mat1);
field_names = fieldnames(params_mat1);

to_appendi = find(field_lengths == n_trs1);

for field_n = 1:length(to_appendi)
    fieldi = to_appendi(field_n);
    curr_field = field_names(fieldi);
    curr_field = curr_field{1};
    
    vec1 = params_mat1.(curr_field);
    vec2 = params_mat2.(curr_field);
    
    params_mat_fused.(curr_field) = [vec1; vec2];
end