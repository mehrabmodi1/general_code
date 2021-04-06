clear all
close all

path = 'C:\Users\Mehrab\Google Drive\Backup\Stuff\Janelia\Talks\DIG\20201009\MB099C_data.xls';  %MBONG2A'1 dataset
gen_path = 'C:\Users\Mehrab\Google Drive\Backup\Stuff\Janelia\Talks\Labmeeting talks\20201218\MB099C_generalization_data.xls'; %Aditya's generalization dataset    
scr_path = 'C:\Users\Mehrab\Google Drive\Backup\Stuff\Janelia\Talks\Labmeeting talks\20201218\Yoshi_cpt_screen.xls';    %Yoshi's US-subst screen across cpts
MB296B_genpath = 'C:\Users\Mehrab\Google Drive\Backup\Stuff\Janelia\Talks\Labmeeting talks\20201218\MB296B_generalization_data.xls';
TH_rescue_discrpath = 'C:\Users\Mehrab\Google Drive\Backup\Stuff\Janelia\Talks\Labmeeting talks\20201218\THRescue_fine_discr.xls';

%Plotting fine discr and coarse discr for multiple compartments
[scr_data] = xlsread(scr_path, 1).*-1;
marker_colors = repmat([0.65, 0.65, 0.65], 10, 1);
line_colors = [];
col_pairs = [];
%Note: this data is for these drivers: MB320C(G1pedc), MB320C, MB099C(G2A'1), MB099C, MB630B(A3), MB630B, MB043C(A1), MB043C, MB213B(B1B2), MB213B,  
xlabels = [{'PA-EL (G1pedc)'}, {'PA-BA'}, {'PA-EL (G2A''1)'}, {'PA-BA'}, {'PA-EL (A3)'}, {'PA-BA'}, {'PA-EL (A1)'}, {'PA-BA'}, {'PA-EL (B1B2)'}, {'PA-BA'}] ;
figure(1)
fig_h = scattered_dot_plot_ttest(scr_data, 1, 1, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 1, [0.8, 0.4, 0.4], 2, 0.05);
hold on
ax_vals = axis;
plot([0, ax_vals(2)], [0, 0], 'Color', [0.65, 0.65, 0.65]);
ylabel('performance index');
fig_wrapup(fig_h, []);


%Plotting fine discr and generalization for MB099C data (G2A'1 + others)
[beh_data] = xlsread(path, 1).*-1;
[gen_data] = xlsread(gen_path, 1);

beh_data = pad_n_concatenate(beh_data, gen_data, 2, nan);
marker_colors = repmat([0.65, 0.65, 0.65], 3, 1);
line_colors = [];
col_pairs = [];
xlabels = [{'PA-EL'}, {'PA-BA'}, {'BA-EL'}];
figure(2)
fig_h = scattered_dot_plot_ttest(beh_data, 2, 1, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 1, [0.8, 0.4, 0.4], 2, 0.05);
ylabel('performance index');
fig_wrapup(fig_h, []);

ax_vals = axis;
ax_vals(4) = 1.1;
axis(ax_vals);
[p, h] = ranksum(beh_data(:, 1), beh_data(:, 2))

[p, h] = ranksum(beh_data(:, 2), beh_data(:, 3))


%Plotting fine discr and generalization for MB296B data (only G2A'1)
[gen_data] = xlsread(MB296B_genpath, 1);
marker_colors = repmat([0.65, 0.65, 0.65], 3, 1);
line_colors = [];
col_pairs = [];
xlabels = [{'PA-EL'}, {'PA-BA'}, {'BA-EL'}];
figure(3)
fig_h = scattered_dot_plot_ttest(gen_data, 3, 1, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 1, [0.8, 0.4, 0.4], 2, 0.05);
ylabel('performance index');
fig_wrapup(fig_h, []);

ax_vals = axis;
ax_vals(4) = 1.1;
axis(ax_vals);


%Plotting TH-rescue, fine discr for MB296B data (only G2A'1)
[rescue_data] = xlsread(TH_rescue_discrpath, 1);
marker_colors = repmat([0.65, 0.65, 0.65], 3, 1);
line_colors = [];
col_pairs = [];
xlabels = [{'G2Ap1rescue'}, {'dbl knckout'}, {'heterozygote'}];
figure(4)
fig_h = scattered_dot_plot_ttest(rescue_data, 4, 1, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 1, [0.8, 0.4, 0.4], 2, 0.05);
ylabel('performance index');
fig_wrapup(fig_h, []);

[p, h] = ranksum(rescue_data(:, 1), rescue_data(:, 2))


ax_vals = axis;
ax_vals(4) = 1.1;
axis(ax_vals);

keyboard
%one-sample ttesting
[p] = signrank_batch(beh_data)    %checking if different from 0
[p] = signrank_batch(scr_data) 
[p] = signrank_batch(rescue_data) 
[p] = signrank_batch(gen_data) 