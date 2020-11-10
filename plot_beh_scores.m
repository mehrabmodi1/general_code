clear all
close all

path = 'C:\Users\Mehrab\Google Drive\Backup\Stuff\Janelia\Talks\DIG\20201009\MB099C_data.xls';  %MBONG2A'1 dataset

[beh_data] = xlsread(path, 1).*-1;

marker_colors = repmat([0.65, 0.65, 0.65], 2, 1);
line_colors = [];
col_pairs = [];
xlabels = [{'PA-EL'}, {'PA-BA'}];
figure(1)
fig_h = scattered_dot_plot_ttest(beh_data, 1, 1, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 1, [0.8, 0.4, 0.4], 1, 0.05);
ylabel('performance index');
fig_wrapup(fig_h, []);

ax_vals = axis;
ax_vals(4) = 1.1;
axis(ax_vals);