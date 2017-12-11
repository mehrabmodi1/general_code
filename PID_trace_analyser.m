clear all
close all

direc = 'D:\Data\Janelia\resonant\20171206\olf_testing2';

param_mat = load([direc, '\params.mat']);
param_mat = param_mat.params_mat;

[stim_mat_simple, column_heads, rand_trains] = load_params_res2(param_mat);
n_trials = size(stim_mat_simple, 1);
odor_list = sort(unique(stim_mat_simple));
n_odors = length(odor_list);
del = find(stim_mat_simple(:, 1) == odor_n);
n_reps = length(del);

for odor_n = 1:n_odors
    curr_trs = find(stim_mat_simple(:, 1) == odor_n);
    
    
    %loading in traces for current odor
    for tr_n = 1:length(curr_trs)
        tr_ni = curr_trs(tr_n);
        trace = load([direc, '\PID_trace_tr-', int2str(tr_ni)]);
        trace = trace.PID_data;

        
    end
    
    
    keyboard
end