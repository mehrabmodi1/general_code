function [stim_mat_simple, column_heads] = generate_stim_mat_simple(stim_mat)
%This function takes an experiment's stimulus specification structure and
%unrolls it into a large, 2-D matrix that is a lot easier to query.

column_heads = [{'odor_n'}, {'duration'}, {'pulse_type'}, {'n_odor_pulses'}, {'inter_pulse_interval'}, {'rand_trains'}, {'mean_rand_pulse_dur'}, {'pulse_train_n'}, {'odour_olf2'}, {'duration_olf2'}, {'rel_stim_latency_olf2'}, {'pulse_type_olf2'}, {'n_od_pulses_olf2'}, {'inter_pulse_interval_olf2'}, {'rand_trains_olf2'}, {'mean_rand_pulse_dur_olf2'}, {'pulse_train_n_olf2'}, {'led_on'}, {'elec_on'}, {'isi'}, {'stim_latency'}, {'post_od_scan_dur'}, {'first_dilution'}, {'second_dilution'}, {'matched_tif_n'}];
n_trials = size(stim_mat, 2);

%loop to read param values from params file into stim_mat
stim_mat_simple = zeros(n_trials, 25) + nan;
for trial_n = 1:n_trials
    stim_mat_simple(trial_n, :) = [stim_mat(trial_n).odours, stim_mat(trial_n).duration, stim_mat(trial_n).pulse_type, stim_mat(trial_n).n_od_pulses, stim_mat(trial_n).inter_pulse_interval, stim_mat(trial_n).rand_trains, stim_mat(trial_n).mean_rand_pulse_dur, stim_mat(trial_n).pulse_train_n, ...
    stim_mat(trial_n).odours_olf2, stim_mat(trial_n).duration_olf2, stim_mat(trial_n).rel_stimLatency_olf2, stim_mat(trial_n).pulse_type_olf2, stim_mat(trial_n).n_od_pulses_olf2, stim_mat(trial_n).inter_pulse_interval_olf2, ...
    stim_mat(trial_n).rand_trains_olf2, stim_mat(trial_n).mean_rand_pulse_dur_olf2,  stim_mat(trial_n).pulse_train_n_olf2, stim_mat(trial_n).led_on, stim_mat(trial_n).elec_on, stim_mat(trial_n).isi, stim_mat(trial_n).stimLatency, stim_mat(trial_n).post_od_scan_dur, stim_mat(trial_n).firstDilution, stim_mat(trial_n).secondDilution, nan];
end
