function [PID_traces] = get_PID_traces(direc, curr_trs, frame_time)
%Syntax:[PID_traces] = get_PID_traces(direc, curr_trs, frame_time)
%This function loads in PID traces, baseline subtracts and normalises them.
%It also down-samples them to match the frame rate.

%testing variables
% direc = 'C:\Data\Data\Analysed_data\20171222\13F02_opGC6f_fluctuating_stim\1\';
% curr_trs = [10 13, 21, 23, 26];
% frame_time = 99.933;        %in ms

n_trials = length(curr_trs);
PID_traces = [];

for trial_n = 1:n_trials
    trial_ni = curr_trs(trial_n);
    curr_trace = load([direc, 'PID_trace_tr-', num2str(trial_ni)]);
    curr_trace = curr_trace.PID_data;
    
    if trial_n == 1
        acqn_time = mean(diff(curr_trace(:, 2)));
        fr_pt_ratio = round(frame_time./acqn_time);
    else
    end
    
    %subtracting baseline
    baseline = mean(curr_trace(1:500, 1));  %this assumes no odor was delivered in the first half second of the trace
    curr_trace(:, 1) = curr_trace(:, 1) - baseline;
    
    %downsampling PID trace to match frame rate
    points_to_sample = linspace(1, max(curr_trace(:, 2)), size(curr_trace, 1)./fr_pt_ratio);
    curr_trace_ds = interp1(curr_trace(:, 2), curr_trace(:, 1), points_to_sample, 'linear');
    
    PID_traces = [PID_traces; curr_trace_ds];
    
end
%normalising while retaining relative differences between repeats
PID_traces = PID_traces./max(max(PID_traces));