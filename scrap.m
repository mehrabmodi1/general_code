clear all
close all

path = 'C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Talks\India_visits\201701\sus_act movies\vlobe_sus1.xls';

data = xlsread(path);
frame_time = .2099;     %in s
stim_time = 16;         %in s
stim_duration = 60;     %in s



stim_frames = [(stim_time./frame_time), (stim_time + stim_duration)./frame_time];

baseline = nanmean(data(1:stim_frames(1)));

dff_data = (data - repmat(baseline, 440, 1))./baseline;


figure(1)
plot(dff_data, 'LineWidth', 3)
add_stim_shading(1, stim_frames, .2, [0.5, 0.7, 0.3]);
set_xlabels_time(1, frame_time, .5)
ylabel('dF/F')