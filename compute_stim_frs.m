function [stim_frs] = compute_stim_frs(stim_mat, trial_n, frame_time)

curr_train = stim_mat(trial_n).rand_trains;      %this applies even for trials that are not random train trials.
stim_latency = stim_mat(trial_n).stim_latency;
stim_frs = compute_pulse_frames_train(curr_train, frame_time, stim_latency);
