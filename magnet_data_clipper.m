function [trace_mat, US_trace_mat, CS_onset_time, CS_duration, CS_US_delay, baseline_scores] = magnet_data_clipper(direc, set_type, control_direc_counter)
% syntax: [trace_mat, CS_onset_times] = magnet_data_clipper(direc, CS_onset_times)
% magnet_data_clipper loads in all initial variables and a matrix of
% magnet traces, appropriately clipped. It takes rand datasets into account
% and yeilds traces assuming that the CS is the stimulus of interest.

%direc = 'C:\Data\data\BlinkView\20120120\field5_rand_assembled\';   %for testing

pre_clip = 2400;     %in ms
post_clip = 2400;    %in ms

pre_clip = pre_clip.*10;    %in no of data points
post_clip = post_clip.* 10; %in no of data points

frame_time = 0.1;           %10 KHz in ms

%initialising dataset specific variables
no_trials = dir([direc '*.tif']);
no_trials = size(no_trials);

tif_tag = imfinfo([direc 'Trial - 0-1.tif']);
img_desc = tif_tag.ImageDescription;

CS_oni = findstr(img_desc, 'CS Onset');
CS_duri = findstr(img_desc, 'CS Dur'); 
CS_US_delayi = findstr(img_desc, 'CS-US');
US_durationi = findstr(img_desc, 'US Dur');
Pb_Tr_i = findstr(img_desc, 'Pb Tr');
CS_minus_i = findstr(img_desc, 'CS- Tr');

%making allowance for rand datasets
if set_type == 1
    CS_onset_times = zeros(1, no_trials) + str2num(img_desc( (CS_oni+16):(CS_oni+19) ));                   %in ms
    US_onset_time = CS_onset_times(1, 1);                                    %CS_US_delay is added below
elseif set_type == 2 && control_direc_counter == 2
    CS_onset_times = load([direc '\rand CS timings.xls']);
    US_onset_time = 12500;
elseif set_type == 2 && control_direc_counter == 1  %for control the datasets in a rand set
    CS_onset_times = zeros(1, no_trials) + str2num(img_desc( (CS_oni+16):(CS_oni+19) ));                   %in ms
    US_onset_time = CS_onset_times(1, 1);
end

CS_duration = str2num(img_desc( (CS_duri+14):(CS_duri+17) ));                   %in ms
CS_US_delay = str2num(img_desc( (CS_US_delayi+20):(CS_US_delayi+23) ));         %in ms
US_onset_time = US_onset_time + CS_US_delay;                                    %in ms
US_duration = str2num(img_desc( (US_durationi+14):(US_durationi+16) ));         %in ms
isProbe_bool = str2num(img_desc( (Pb_Tr_i+16)) );                               %logical
isCSminus_bool = str2num(img_desc( (CS_minus_i+17)) );                          %logical


US_onset_frame = US_onset_time./frame_time;


trace = load([direc '\Trial - 0-1.xls']);
no_frames = size(trace, 1);


clear trace
clear CS_oni
clear CS_duri
clear CS_US_delayi
clear US_durationi
clear PB_TR_i
clear CS_minus_i

%loading in traces
trace_mat = zeros(( (pre_clip + 6000 + post_clip)), no_trials);
US_trace_mat = zeros( (600./frame_time), no_trials);
baseline_scores = zeros(no_trials, 1);
for trial_no = 1:no_trials
    CS_onset_time = CS_onset_times(1, trial_no);                %in ms
    CS_onset_frame = CS_onset_time./frame_time;
    trace = load([direc 'Trial - ' int2str(trial_no-1) '-1.xls']);
    trace = trace(:, 1);
    if CS_onset_frame > pre_clip && (CS_onset_frame + ( (600./frame_time + post_clip)) - 1) < size(trace, 1)          %in case rand CS_onset_time too close to either end of trace
        baseline_trace = trace( (CS_onset_frame - pre_clip):(CS_onset_frame - 1), 1 );
        baseline_score = median(baseline_trace);
        US_baseline = trace((US_onset_frame - (600./frame_time)):(US_onset_frame - 1) );
        US_baseline = median(US_baseline);
        
        US_trace = trace(US_onset_frame:(US_onset_frame + (600./frame_time) - 1), 1) - US_baseline;
        trace = trace((CS_onset_frame - (pre_clip) ):(CS_onset_frame + ( (600./frame_time + post_clip)) - 1), 1) - baseline_score;
    elseif CS_onset_frame < pre_clip
        baseline_trace = trace(1:(CS_onset_frame - 1), 1 );
        baseline_score = median(baseline_trace);
        pad = zeros( ((pre_clip - CS_onset_frame) + 1), 1);
        US_baseline = trace((US_onset_frame - (600./frame_time)):(US_onset_frame - 1) );
        US_baseline = median(US_baseline);
        
        US_trace = trace(US_onset_frame:(US_onset_frame + (600./frame_time) - 1), 1) - US_baseline;
        trace = [pad; trace];
        trace = trace(( (CS_onset_frame + length(pad)) - (pre_clip) ):( (CS_onset_frame + length(pad)) + ( (600./frame_time + post_clip)) - 1), 1) - baseline_score;
    elseif (CS_onset_frame + ( (600./frame_time + post_clip)) - 1) > size(trace, 1)
        baseline_trace = trace( (CS_onset_frame - pre_clip):(CS_onset_frame - 1), 1 );
        baseline_score = median(baseline_trace);
        US_baseline = trace((US_onset_frame - (600./frame_time)):(US_onset_frame - 1) );
        US_baseline = median(US_baseline);
        
        US_trace = trace(US_onset_frame:(US_onset_frame + (600./frame_time) - 1), 1) - US_baseline;
        pad = zeros((CS_onset_frame + (600./frame_time) + post_clip) - size(trace, 1), 1);
        trace = [trace; pad];
        trace = trace((CS_onset_frame - (pre_clip) ):(CS_onset_frame + ( (600./frame_time + post_clip) ) - 1), 1) - baseline_score;
    end
    trace_mat(:, trial_no) = trace;
    US_trace_mat(:, trial_no) = US_trace;
    baseline_scores(trial_no, 1) = baseline_score;
    CS_onset_time = pre_clip.*frame_time;
end