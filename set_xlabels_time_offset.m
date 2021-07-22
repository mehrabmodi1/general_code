function [] = set_xlabels_time_offset(fig_n, frame_time, t_step_multiplier, t_offset)
%syntax:[] = set_xlabels_time(fig_n, frame_time, t_step_multiplier, t_offset)
%Note: Specify frame_time in s.

figure(fig_n)
axis_m = axis;

n_frames = floor(axis_m(1, 2));
tot_time = n_frames.*frame_time;
range_n = floor(log10(tot_time));   %10 raised to range_n = tot_time, useful for rounding off t_step
t_step = 10.^(range_n - 1).*t_step_multiplier;

zero_val = 0 + t_offset;
tot_time_orig = tot_time;
tot_time = tot_time + t_offset;

x_tick_labelsi = zero_val:t_step:tot_time; 
x_ticks = 0:t_step:tot_time_orig./frame_time;

ax = gca;
ax.XTick = x_ticks;

%building list of time point strings
for tick_n = 1:length(x_tick_labelsi)
    x_tick_labels{1, tick_n} = num2str(x_tick_labelsi(tick_n));
end
ax.XTickLabel = x_tick_labels;
xlabel('time (s)')
