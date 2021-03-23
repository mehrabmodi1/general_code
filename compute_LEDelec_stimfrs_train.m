function [stim_frs_train, stim_frs_train_rounded] = compute_LEDelec_stimfrs_train(stim_mat, tr_n, frame_time)
%Syntax:[stim_frs_train, stim_frs_train_rounded] = compute_LEDelec_stimfrs_train(stim_mat, tr_n, frame_time)
%This function computes the frame numbers when each electrical or LED
%stimulus pulse was delivered on a specified trial. Since these pulses are
%usually short ( < 1 frame time), a rounded off version of the pulse train
%is also computed.

stim_onset_t = stim_mat(tr_n).rel_stim_init_delay + stim_mat(tr_n).stimLatency;     %train onset time in s
pulseon_width = (1./stim_mat(tr_n).stim_freq).*(stim_mat(tr_n).st_duty_cyc./100);   %in s
pulseoff_width = (1./stim_mat(tr_n).stim_freq).*(1 - (stim_mat(tr_n).st_duty_cyc./100));   %in s
n_pulses = stim_mat(tr_n).stim_dur.*stim_mat(tr_n).stim_freq;

stim_frs_train = [];
next_pulse_t = stim_onset_t;
for pulse_n = 1:n_pulses
    stim_frs_train = [stim_frs_train; [next_pulse_t, (next_pulse_t + pulseon_width)]];
    next_pulse_t = next_pulse_t + pulseon_width + pulseoff_width;
end

stim_frs_train = stim_frs_train./frame_time;
stim_frs_train_rounded = [floor(stim_frs_train(:, 1)), ceil(stim_frs_train(:, 2))];

end