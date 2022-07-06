clear all
close all

n_pulses = 10;
% on_durs = abs(randn(n_pulses, 1));
% off_durs = abs(randn(n_pulses, 1)).*2;
on_durs = 2.^(0:0.5:5)';
off_durs = max(on_durs).*1.2 - on_durs;

off_durs = [2; off_durs];

pulse_times = [0; 2];
pulse_levels = [0; 0];
for pulse_n = 1:n_pulses
    
    off_pt = pulse_times(end) + off_durs(pulse_n);
    pulse_times = [pulse_times; off_pt; off_pt]
    pulse_levels = [pulse_levels; 0; 1];
    
    on_pt = off_pt + on_durs(pulse_n);
    pulse_times = [pulse_times; on_pt; on_pt]
    pulse_levels = [pulse_levels; 1; 0];    
    
end

plot(pulse_times, pulse_levels);
ax_vals = axis;
ax_vals(4) = 1.2;
axis(ax_vals);


