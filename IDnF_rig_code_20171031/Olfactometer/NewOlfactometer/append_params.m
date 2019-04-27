function [params_struc_fused] = append_params(params_mat1, params_mat2, randomise)

n_trs1 = size(params_mat1, 2);
n_trs2 = size(params_mat2, 2);


params_struc_fused = [params_mat1, params_mat2];

%randomising trials in fused param matrix if specified
if randomise == 1
    params_struc_fused_orig = params_struc_fused;
    n_trs = n_trs1 + n_trs2;
    r_vec = randperm(n_trs);
    for tr_n = 1:n_trs
        params_struc_fused(tr_n) = params_struc_fused_orig(r_vec(tr_n));
    end 
        
else
end
