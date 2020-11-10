function [stim_frs, bigger_stim_olf] = compute_stim_frs_modular(stim_mat, trial_n, frame_time)
%syntax: [stim_frs, bigger_stim_olf] = compute_stim_frs_modular(stim_mat, trial_n, frame_time)

%computing stim_frs for olfactometer1
curr_train = stim_mat(trial_n).pulse_train;      %this applies even for trials that are not random train trials.
stim_latency = stim_mat(trial_n).stimLatency;

stim_frs_olf1 = compute_pulse_frames_train(curr_train, frame_time, stim_latency);

if isnan(stim_frs_olf1) == 0
    tot_stim_area_olf1 = sum(stim_frs_olf1(:, 2) - stim_frs_olf1(:, 1));
else
    tot_stim_area_olf1 = 0;     %case when olfactometer2 was not used
end
%computing stim_frs for olfactometer2
curr_train = stim_mat(trial_n).pulse_train_olf2;      %this applies even for trials that are not random train trials.
try
    stim_latency = stim_mat(trial_n).rel_stimLatency_olf2 + stim_latency;
catch
    keyboard
end
stim_frs_olf2 = compute_pulse_frames_train(curr_train, frame_time, stim_latency);
if isnan(stim_frs_olf2) == 0
    tot_stim_area_olf2 = sum(stim_frs_olf2(:, 2) - stim_frs_olf2(:, 1));
else
    tot_stim_area_olf2 = 0;     %case when olfactometer2 was not used
end

stim_frs = [{stim_frs_olf1}, {stim_frs_olf2}];

[del, bigger_stim_olf] = max([tot_stim_area_olf1, tot_stim_area_olf2]);     %identifying olfactometer number with more stimulus area
