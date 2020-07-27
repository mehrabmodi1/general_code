function [PID_traces, traces_orig] = get_PID_traces(direc, curr_trs, frame_time, normalise)
%Syntax:[PID_traces] = get_PID_traces(direc, curr_trs, frame_time)
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
      
    if trial_n == 1
        acqn_time = mean(diff(curr_trace(:, 3)));
        fr_pt_ratio = round(frame_time./acqn_time);
    else
    end
    
    %subtracting baseline
    baseline = mean(curr_trace(400:500, 1:2), 1);  %this assumes no odor was delivered in the first half second of the trace
    curr_trace(:, 1:2) = curr_trace(:, 1:2) - repmat(baseline, size(curr_trace, 1), 1);
    
    %downsampling PID trace to match frame rate
    x = curr_trace(:, 1);
    n = fr_pt_ratio;
    m = numel(x);
    curr_trace_ds(:, 1) = mean(reshape( [x(:);nan(mod(-m,n),1)],n,[]), 'omitnan');
   
    x = curr_trace(:, 2);
    n = fr_pt_ratio;
    m = numel(x);
    curr_trace_ds(:, 2) = mean(reshape( [x(:);nan(mod(-m,n),1)],n,[]), 'omitnan');
    curr_trace_ds(1:4, :) = 0;
    %curr_trace_ds = interp1(curr_trace(:, 3), curr_trace(:, 1), points_to_sample, 'linear');
    PID_traces = [PID_traces, curr_trace_ds];
    try
        traces_orig = [traces_orig, curr_trace];
    catch
        keyboard
    end
end

if normalise == 1
    %normalising while retaining relative differences between repeats
    PID_traces = PID_traces./max(max(PID_traces));
else
end