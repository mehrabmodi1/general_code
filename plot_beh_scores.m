clear all
close all


scr_path = 'C:\Data\Data\beh_data_extr_PIs\Yoshi_cpt_screen.xls';    %Yoshi's US-subst screen across cpts
A3_path = 'C:\Data\Data\beh_data_extr_PIs\Alpha3_fine_coarse.xls';   %US-subst 3x training, 1 day memory 
MB099C_fine_path = 'C:\Data\Data\beh_data_extr_PIs\MB099C_fine_discr.xls';  %MBONG2A'1 dataset
MB099C_gen_path = 'C:\Data\Data\beh_data_extr_PIs\MB099C_gen_data.xls'; %Aditya's generalization dataset    
MB296B_genpath = 'C:\Data\Data\beh_data_extr_PIs\MB296B_gen_data.xls';
TH_rescue_discrpath = 'C:\Users\Mehrab\Google Drive\Backup\Stuff\Janelia\Talks\Labmeeting talks\20201218\THRescue_fine_discr.xls';

%Plotting fine discr and coarse discr for multiple compartments
[scr_data] = xlsread(scr_path, 1).*-1;
marker_colors = repmat([0.65, 0.65, 0.65], 10, 1);
line_colors = [];
col_pairs = [];
%Note: this data is for these drivers: MB320C(G1pedc), MB320C, MB099C(G2A'1), MB099C, MB043C(A1), MB043C, MB213B(B1B2), MB213B,  
xlabels = [{'PA-EL (G1pedc)'}, {'PA-BA'}, {'PA-EL (G2A''1)'}, {'PA-BA'}, {'PA-EL (A1)'}, {'PA-BA'}, {'PA-EL (B1B2)'}, {'PA-BA'}] ;

%normalizing to easy discr median PIs
% median_vec = median(scr_data(:, [1, 3, 5, 7]), 'omitnan');
% median_vec = interleave2(median_vec, median_vec, 'col');
% scr_data = scr_data./repmat(median_vec, size(scr_data, 1), 1);

figure(1)
fig_h = scattered_dot_plot_ttest(scr_data, 1, 1, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 2, [0.8, 0.4, 0.4], 2, 0.05);
hold on
ax_vals = axis;
plot([0, ax_vals(2)], [0, 0], 'Color', [0.65, 0.65, 0.65]);
ylabel('performance index');
fig_wrapup(fig_h, []);


%Plotting fine, coarse discr data for Alpha3, 3x training, 1 day memory
[alpha3_data] = xlsread(A3_path, 1).*-1;
marker_colors = repmat([0.65, 0.65, 0.65], 2, 1);
line_colors = [];
col_pairs = [];
%Note: this data is for these drivers: MB630B(B1B2), MB630B,  
xlabels = [{'PA-EL (A3)'}, {'PA-BA'}] ;

figure(2)
fig_h = scattered_dot_plot_ttest(alpha3_data, 2, 1, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 2, [0.8, 0.4, 0.4], 2, 0.05);
hold on
ax_vals = axis;
plot([0, ax_vals(2)], [0, 0], 'Color', [0.65, 0.65, 0.65]);
ylabel('1 day memory PI');
fig_wrapup(fig_h, []);

%Plotting fine discr and generalization for MB099C data (G2A'1 + others)
[beh_data] = xlsread(MB099C_fine_path, 1).*-1;
beh_data_MB099C = beh_data;
[gen_data] = xlsread(MB099C_gen_path, 1);

beh_data = pad_n_concatenate(beh_data, gen_data, 2, nan);
left_color = [0.35, 0.75, 0.35];
right_color = [0.45, 0.45, 0.85];
marker_colors = repmat(left_color, 2, 1);
marker_colors = [marker_colors; repmat(right_color, 2, 1)];

line_colors = [];
col_pairs = [];
xlabels = [{'PA-EL'}, {'PA-BA'}, {'PA-EL'}, {'BA-EL'}];
figure(3)
fig_h = scattered_dot_plot_ttest(beh_data, 3, 1, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 2, [0.8, 0.4, 0.4], 2, 0.05);
ylabel('discrimination score');
ax_vals = axis;
ax_vals(4) = 1.1;
axis(ax_vals);
ax_vals = axis;
yyaxis('right');
axis(ax_vals);
ylabel('generalization score');
ax = gca;
ax.YAxis(1).Color = left_color;
ax.YAxis(2).Color = right_color;
fig_wrapup(fig_h, []);


[p, h] = ranksum(beh_data(:, 1), beh_data(:, 2))

[p, h] = ranksum(beh_data(:, 2), beh_data(:, 3))


%Plotting fine discr and generalization for MB296B data (only G2A'1)
beh_data = beh_data(:, 1:2).*0.0001;
[gen_data_MB296B] = xlsread(MB296B_genpath, 1);
MB296B_data = pad_n_concatenate(beh_data, gen_data_MB296B, 2, nan);
left_color = [0.35, 0.75, 0.35];
right_color = [0.45, 0.45, 0.85];
marker_colors = repmat(left_color, 2, 1);
marker_colors = [marker_colors; repmat(right_color, 2, 1)];
line_colors = [];
col_pairs = [];
xlabels = [{'PA-EL'}, {'PA-BA'}, {'PA-EL'}, {'BA-EL'}];
figure(4)
fig_h = scattered_dot_plot_ttest(MB296B_data, 4, 1, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 2, [0.8, 0.4, 0.4], 2, 0.05);
ylabel('discrimination score');
ax_vals = axis;
ax_vals(4) = 1.1;
axis(ax_vals);
ax_vals = axis;
yyaxis('right');
axis(ax_vals);
ylabel('generalization score');
ax = gca;
ax.YAxis(1).Color = left_color;
ax.YAxis(2).Color = right_color;
fig_wrapup(fig_h, []);

ax_vals = axis;
ax_vals(4) = 1.1;
axis(ax_vals);


%Plotting TH-rescue, fine discr for MB296B data (only G2A'1)
[rescue_data] = xlsread(TH_rescue_discrpath, 1);
marker_colors = repmat([0.65, 0.65, 0.65], 2, 1);
line_colors = [];
col_pairs = [];
xlabels = [{'G2Ap1rescue'}, {'dbl knckout'}];
figure(5)
fig_h = scattered_dot_plot_ttest(rescue_data(:, [1:2]), 5, 1, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 2, [0.8, 0.4, 0.4], 2, 0.05);
ylabel('performance index');
fig_wrapup(fig_h, []);

[p, h] = ranksum(rescue_data(:, 1), rescue_data(:, 2))


ax_vals = axis;
ax_vals(4) = 1.1;
axis(ax_vals);

%one-sample ttesting
[p_MB099C_discr] = signrank_batch(beh_data)    %checking if different from 0
[p_MB099C_gen] = signrank_batch(gen_data)
[p_screen] = signrank_batch(scr_data) 
[p_rescue] = signrank_batch(rescue_data(:, [1:2])) 
 