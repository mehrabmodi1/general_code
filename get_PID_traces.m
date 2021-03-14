function [PID_traces, traces_orig, n_data_dims] = get_PID_traces(direc, curr_trs, frame_time, normalise)
%Syntax:[PID_traces] = get_PID_traces(direc, curr_trs, frame_time, normalise)
%This function loads in PID traces, baseline subtracts and normalises them.
%It also down-samples them to match the frame rate.

%testing variables
% direc = 'C:\Data\Data\Analysed_data\20171222\13F02_opGC6f_fluctuating_stim\1\';
% curr_trs = [10 13, 21, 23, 26];
% frame_time = 99.933;        %in ms

%frame_time = frame_time.*1000;  %converting from s to ms
n_trials = length(curr_trs);
PID_traces = [];
traces_orig = [];
for trial_n = 1:n_trials
    trial_ni = curr_trs(trial_n);
    curr_trace = load([direc, 'PID_trace_tr-', num2str(trial_ni)]);
    curr_trace = curr_trace.PID_data;
    
    n_data_dims = size(curr_trace, 2) - 1;
    
    if trial_n == 1
        acqn_time = mean(diff(curr_trace(:, size(curr_trace, 2))));       %computing Analog acqn time; column 3 is time of each sample in this case
        fr_pt_ratio = round(frame_time./acqn_time);
    else
    end
  
    %subtracting baseline
    baseline = mean(curr_trace(400:500, 1:n_data_dims), 1);  %this assumes no odor was delivered in the first half second of the trace
    curr_trace(:, 1:n_data_dims) = curr_trace(:, 1:n_data_dims) - repmat(baseline, size(curr_trace, 1), 1);
    
    %downsampling PID trace to match frame rate
    n = fr_pt_ratio;
    curr_trace_ds = [];
    
    for data_dim_n = 1:n_data_dims
        x = curr_trace(:, data_dim_n);
        m = numel(x);
        curr_trace_ds(:, data_dim_n) = mean(reshape( [x(:);nan(mod(-m,n),1)],n,[]), 'omitnan')';
    end
    
    curr_trace_ds(1:4, :) = 0;
    %curr_trace_ds = interp1(curr_trace(:, 3), curr_trace(:, 1), points_to_sample, 'linear');
    PID_traces = pad_n_concatenate(PID_traces, curr_trace_ds, 2, nan);
    try
        traces_orig = pad_n_concatenate(traces_orig, curr_trace, 2, nan);
    catch
        keyboard
    end
end

if normalise == 1
    %normalising while retaining relative differences between repeats
    PID_traces = PID_traces./max(max(PID_traces));
else
end