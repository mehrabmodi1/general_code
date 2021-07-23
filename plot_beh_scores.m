clear all
close all


scr_path = 'C:\Data\Data\beh_data_extr_PIs\Yoshi_cpt_screen.xls';    %Yoshi's US-subst screen across cpts
A3_path = 'C:\Data\Data\beh_data_extr_PIs\Alpha3_fine_coarse.xls';   %US-subst 10x training, 1 day memory 
A3_genpath = 'C:\Data\Data\beh_data_extr_PIs\Alpha3_generalization.xls';   %US-subst 10x training, 1 day memory 
MB099C_fine_path = 'C:\Data\Data\beh_data_extr_PIs\MB099C_fine_discr.xls';  %MBONG2A'1 dataset
MB099C_gen_path = 'C:\Data\Data\beh_data_extr_PIs\MB099C_gen_data.xls'; %Aditya's generalization dataset    
MB296B_genpath = 'C:\Data\Data\beh_data_extr_PIs\MB296B_gen_data_3xtraining.xls';
MB296B_discrpath = 'C:\Data\Data\beh_data_extr_PIs\MB296B_discr_data_3xtraining.xls';
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

mean_color = [0, 0, 0];

wrapup_vars{1} = [100, 25];
wrapup_vars{2} = 0.6;

figure(1)
fig_h = scattered_dot_plot_ttest(scr_data, 1, .6, 1, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 2, 0.05, 0, 1, 'force_mean', wrapup_vars);
hold on
ax_vals = axis;
ax_vals(3) = -0.8;
ax_vals(4) = 0.8;
axis(ax_vals);
plot([0, ax_vals(2)], [0, 0], 'Color', [0, 0, 0]);
ylabel('immediate discrimination (PI)');
fig_wrapup(fig_h, []);

%plotting only aversive cpts
xlabels = [{'PA-EL (G1)'}, {'PA-BA'}];
figure(6)
fig_h = scattered_dot_plot_ttest(scr_data(:, 1:2), 6, .6, 1, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 2, 0.05, 0, 1, 'force_mean', wrapup_vars);
hold on
ax_vals = axis;
ax_vals(3) = -0.2;
ax_vals(4) = 0.85;
axis(ax_vals);
plot([0, ax_vals(2)], [0, 0], 'Color', [0, 0, 0]);
ylabel('imm. discrimination (PI)');
fig_wrapup(fig_h, [], [100, 25], 0.6);


%plotting only aversive cpts
xlabels = [{'PA-EL (G2)'}, {'PA-BA'}];
figure(9)
fig_h = scattered_dot_plot_ttest(scr_data(:, 3:4), 9, .6, 1, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 2, 0.05, 0, 1, 'force_mean', wrapup_vars);
hold on
ax_vals = axis;
ax_vals(3) = -0.2;
ax_vals(4) = 0.85;
axis(ax_vals);
plot([0, ax_vals(2)], [0, 0], 'Color', [0, 0, 0]);
ylabel('imm. discrimination (PI)');
fig_wrapup(fig_h, [], [100, 25], 0.6);




%Plotting fine, coarse discr data for Alpha3, 3x training, 1 day memory
[alpha3_data] = xlsread(A3_path, 1).*-1;
marker_colors = repmat([0.65, 0.65, 0.65], 2, 1);
line_colors = [];
col_pairs = [];
%Note: this data is for these drivers: MB630B(B1B2), MB630B,  
xlabels = [{'PA-EL (A3)'}, {'PA-BA'}] ;

figure(2)
fig_h = scattered_dot_plot_ttest(alpha3_data, 2, .6, 1, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 2, 0.05, 0, 1, 'force_mean', wrapup_vars);
hold on
ax_vals = axis;
ax_vals(3) = -0.2;
ax_vals(4) = 0.85;
axis(ax_vals);
plot([0, ax_vals(2)], [0, 0], 'Color', [0. 0. 0]);
ylabel('1 day memory PI');
fig_wrapup(fig_h, [], [100, 25], 0.6);

%Plotting fine discr and generalization for MB099C data (G2A'1 + others)
[beh_data] = xlsread(MB099C_fine_path, 1).*-1;
beh_data_MB099C = beh_data;
[gen_data] = xlsread(MB099C_gen_path, 1);

beh_data = pad_n_concatenate(beh_data, gen_data, 2, nan);
left_color = [0.65, 0.65, 0.65];
right_color = [0, 0, 0];
marker_colors = repmat(left_color, 2, 1);
marker_colors = [marker_colors; repmat(right_color, 2, 1)];

line_colors = [];
col_pairs = [];
xlabels = [{'PA-EL'}, {'PA-BA'}, {'PA-EL'}, {'BA-EL'}];
figure(3)
fig_h = scattered_dot_plot_ttest(beh_data, 3, .6, 1, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 2, 0.05, 0, 1, 'force_mean', wrapup_vars);
ylabel('discrimination score');
ax_vals = axis;
ax_vals(4) = 0.85;
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
[discr_data_MB296B] = xlsread(MB296B_discrpath, 1);
[gen_data_MB296B] = xlsread(MB296B_genpath, 1);
MB296B_data = pad_n_concatenate(discr_data_MB296B, gen_data_MB296B, 2, nan);
% left_color = [0.35, 0.75, 0.35];
% right_color = [0.45, 0.45, 0.85];
marker_colors = repmat(left_color, 2, 1);
marker_colors = [marker_colors; repmat(right_color, 2, 1)];
line_colors = [];
col_pairs = [];
xlabels = [{'PA-EL'}, {'PA-BA'}, {'PA-EL'}, {'BA-EL'}];
figure(4)
fig_h = scattered_dot_plot_ttest(MB296B_data, 4, .6, 1, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 2, 0.05, 0, 1, 'force_mean', wrapup_vars);
ylabel('discrimination score (3x)');
ax_vals = axis;
ax_vals(4) = 0.85;
axis(ax_vals);
ax_vals = axis;
yyaxis('right');
axis(ax_vals);
ylabel('generalization score (3x)');
ax = gca;
ax.YAxis(1).Color = left_color;
ax.YAxis(2).Color = right_color;
fig_wrapup(fig_h, []);

ax_vals = axis;
ax_vals(4) = 0.85;
axis(ax_vals);

%Plotting MB296B discr and gen scores separately
figure(7)
xlabels = [{'G2 PA-EL'}, {'BA-EL'}];
fig_h = scattered_dot_plot_ttest(discr_data_MB296B, 7, .6, 1, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 2, 0.05, 0, 1, 'force_mean', wrapup_vars);
hold on
ylabel('discrimination score (PI)');
ax_vals = axis;
ax_vals(3) = -0.2;
ax_vals(4) = 0.85;
axis(ax_vals);
plot([0, ax_vals(2)], [0, 0], 'Color', [0, 0, 0]);
fig_wrapup(fig_h, [], [100, 25], 0.6);

figure(8)
xlabels = [{'G2PAEL'}, {'BAEL'}];
fig_h = scattered_dot_plot_ttest(gen_data_MB296B, 8, .6, 1, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 2, 0.05, 0, 1, 'force_mean', wrapup_vars);
hold on
plot([0, ax_vals(2)], [0, 0], 'Color', [0, 0, 0]);
ylabel('generalization score');
ax_vals = axis;
ax_vals(4) = 0.85;
axis(ax_vals);
fig_wrapup(fig_h, [], [100, 25], 0.6);

[gen_data_Alpha3] = xlsread(A3_genpath, 1);
figure(9)
xlabels = [{'A3PAEL'}, {'BAEL'}];
fig_h = scattered_dot_plot_ttest(gen_data_Alpha3, 8, .6, 1, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 2, 0.05, 0, 1, 'force_mean', wrapup_vars);
hold on
plot([0, ax_vals(2)], [0, 0], 'Color', [0, 0, 0]);
ylabel('generalization score');
ax_vals = axis;
ax_vals(4) = 0.85;
axis(ax_vals);
fig_wrapup(fig_h, [], [100, 25], 0.6);


%Plotting TH-rescue, fine discr for MB296B data (only G2A'1)
[rescue_data] = xlsread(TH_rescue_discrpath, 1);
marker_colors = repmat([0.65, 0.65, 0.65], 2, 1);
line_colors = [];
col_pairs = [];
xlabels = [{'G2res'}, {'dblknk'}];
figure(5)
fig_h = scattered_dot_plot_ttest(rescue_data(:, [1:2]), 5, .6, 1, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 2, 0.05, 0, 1, 'force_mean', wrapup_vars);
ylabel('performance index');
fig_wrapup(fig_h, [], [100, 25], 0.6);

[p, h] = ranksum(rescue_data(:, 1), rescue_data(:, 2))


ax_vals = axis;
ax_vals(4) = 0.85;
axis(ax_vals);

%one-sample ttesting
[p_MB099C_discr] = signrank_batch(beh_data)    %checking if different from 0
[p_MB099C_gen] = signrank_batch(gen_data)
[p_MB296B_discr] = signrank_batch(discr_data_MB296B)
[p_MB296B_gen] = signrank_batch(gen_data_MB296B)
[p_screen] = signrank_batch(scr_data) 
[p_rescue] = signrank_batch(rescue_data(:, [1:2])) 
 p_aplha3 = signrank_batch(alpha3_data(:, [1:2])) 