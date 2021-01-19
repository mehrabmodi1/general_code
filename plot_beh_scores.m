clear all
close all

path = 'C:\Users\Mehrab\Google Drive\Backup\Stuff\Janelia\Talks\DIG\20201009\MB099C_data.xls';  %MBONG2A'1 dataset
gen_path = 'C:\Users\Mehrab\Google Drive\Backup\Stuff\Janelia\Talks\Labmeeting talks\20201218\MB099C_generalization_data.xls'; %Aditya's generalization dataset    

[beh_data] = xlsread(path, 1).*-1;
[gen_data] = xlsread(gen_path, 1);

beh_data = pad_n_concatenate(beh_data, gen_data, 2, nan);
marker_colors = repmat([0.65, 0.65, 0.65], 3, 1);
line_colors = [];
col_pairs = [];
xlabels = [{'PA-EL'}, {'PA-BA'}, {'BA-EL'}];
figure(1)
fig_h = scattered_dot_plot_ttest(beh_data, 1, 1, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 1, [0.8, 0.4, 0.4], 1, 0.05);
ylabel('performance index');
fig_wrapup(fig_h, []);

ax_vals = axis;
ax_vals(4) = 1.1;
axis(ax_vals);

%one-sample ttesting
[h, p] = ttest(beh_data)