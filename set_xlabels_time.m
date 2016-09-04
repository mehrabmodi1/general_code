function [] = set_xlabels_time(fig_n, frame_time, n_ticks_multiplier)

figure(fig_n)
axis_m = axis;

n_frames = floor(axis_m(1, 2));
tot_time = n_frames.*frame_time;
range_n = floor(log10(tot_time));   %10 raised to range_n = tot_time, useful for rounding off t_step
tot_time_rounded = round(tot_time, (range_n-1));
n_ticks = round(((tot_time_rounded./(10.^range_n)).*2).*n_ticks_multiplier);

t_step = round(tot_time_rounded./n_ticks, (range_n-1) );

x_tick_labelsi = 0:t_step:tot_time;
x_ticks = x_tick_labelsi./frame_time;

ax = gca;
ax.XTick = x_ticks;

%building list of time point strings
for tick_n = 1:length(x_tick_labelsi)
    x_tick_labels{1, tick_n} = num2str(x_tick_labelsi(tick_n));
end
ax.XTickLabel = x_tick_labels;
xlabel('time (s)')