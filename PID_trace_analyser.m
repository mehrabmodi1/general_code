clear all
close all

direc = 'D:\Data\Janelia\resonant\20171213\olf_calib3';

param_mat = load([direc, '\params.mat']);
param_mat = param_mat.params_mat;

[stim_mat_simple, column_heads, rand_trains] = load_params_res2(param_mat);
n_trials = size(stim_mat_simple, 1);
odor_list = sort(unique(stim_mat_simple(:, 1)));
n_odors = length(unique(odor_list));
n_rand_trains = max(stim_mat_simple(:, 10));
del = find(stim_mat_simple(:, 1) == 1);
n_reps = length(del);
stim_point = round(stim_mat_simple(1, 6).*2000);        %PID datapoint number when stim begins. stim_latency multiplied by the acq rate ie. 2KHz
stim_end_point = round(stim_mat_simple(1, 2)).*2000 + stim_point;

mean_trace_mat = [];
for odor_n = 1:n_odors
    odor_ni = odor_list(odor_n);
    
    for rand_train_n = 1:n_rand_trains
        curr_trs = find(stim_mat_simple(:, 1) == odor_ni & stim_mat_simple(:, 10) == rand_train_n);
        trace_mat = [];
        %loading in traces for current odor
        for tr_n = 1:length(curr_trs)
            tr_ni = curr_trs(tr_n);
            trace = load([direc, '\PID_trace_tr-', int2str(tr_ni)]);
            trace = trace.PID_data;

            if tr_n == 1
                trace_mat = [trace_mat, trace(:, 1)];
                time_vec = trace(:, 2);
            elseif tr_n > 1
                trace_mat = pad_n_concatenate(trace_mat, trace(:, 1), 2, nan);
            end
        end
    
        %calculating dF/F
        mean_baselines = mean(trace_mat(1:stim_point, :), 1, 'omitnan');
        mean_baselines = repmat(mean_baselines, size(trace_mat, 1), 1);
        trace_mat = (trace_mat - mean_baselines);
        
        trace_mat_all(:, :, rand_train_n, odor_n) = trace_mat;
        mean_trace_mat = [mean_trace_mat, mean(trace_mat, 2, 'omitnan')];
    
    end
    

end
plot(time_vec, squeeze(trace_mat_all(:, :, 1, 1)))

%Analysing traces
init_pks = max(mean_trace_mat(stim_point:(stim_point + 10000), :)); 
late_pks = max(mean_trace_mat((stim_end_point - 10000):stim_end_point, :)); 

percentage_drops = (init_pks - late_pks)./init_pks;

