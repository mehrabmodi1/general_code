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
acq_point = round(stim_mat_simple(1, 6).*2000);        %PID datapoint number when stim begins. stim_latency multiplied by the acq rate ie. 2KHz

trace_mat_all = [];
for odor_n = 1:n_odors
    odor_ni = odor_list(odor_n);
    curr_trs = find(stim_mat_simple(:, 1) == odor_ni);
    for rand_train_n = 1:n_rand_trains
        trace_mat = [];
        %loading in traces for current odor
        for tr_n = 1:length(curr_trs)
            tr_ni = curr_trs(tr_n);
            trace = load([direc, '\PID_trace_tr-', int2str(tr_ni)]);
            trace = trace.PID_data;

            if tr_n == 1
                trace_mat = [trace_mat, trace(:, 1)];
            elseif tr_n > 1
                trace_mat = pad_n_concatenate(trace_mat, trace(:, 1), 2, nan);
            end
        end
    
        %calculating dF/F
        mean_baselines = mean(trace_mat(1:acq_point, :), 1, 'omitnan');

        mean_baselines = repmat(mean_baselines, size(trace_mat, 1), 1);
        trace_mat = (trace_mat - mean_baselines);

        if odor_n == 1
            trace_mat_all(:, :, 1) = trace_mat;
            mean_trace_mat = [mean_trace_mat, mean(trace_mat, 2, 'omitnan')];
        elseif odor_n > 1
            trace_mat_all = pad_n_concatenate(trace_mat_all, trace_mat, 3, nan);
        end
    end

end
plot(squeeze(trace_mat_all(:, :, 2)))

%Analysing traces
init_pks = 
