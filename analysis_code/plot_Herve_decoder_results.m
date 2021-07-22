clear all
close all

path_base = 'C:\Users\Mehrab\Google Drive\Backup\Stuff\Janelia\paper_drafts\Mehrab_papers\PaBaEl2\images\figures\figure2_KCs_simple\figures_Herve\round3\extr_data\';


%plotting perceptron output
%1. Gamma, fly4
curr_path = [path_base, 'resps-gamma_4_singles.hdf5'];
curr_metadata  = h5info(curr_path);
curr_data = h5read(curr_path, '/resps');
curr_tr_list = h5read(curr_path, '/odor_tri');

plot_decoder_outputs(curr_data, curr_tr_list);

keyboard

close figure 1
close figure 2
close figure 3

%2. ApBp, fly2
curr_path = [path_base, 'resps-alpha_pbeta_p_2_singles.hdf5'];
curr_metadata  = h5info(curr_path);
curr_data = h5read(curr_path, '/resps');
curr_tr_list = h5read(curr_path, '/odor_tri');

plot_decoder_outputs(curr_data, curr_tr_list);


%--------------------------------

%Plotting perceptron accuracies




function [] = plot_decoder_outputs(curr_data, curr_tr_list)

frame_time = 0.099;
stim_time = [10, 15];   %in s, the od on, od off times
stim_frs = round(stim_time./frame_time);

color_vecs = [0, 0, 0; 0.65, 0.65, 0.65; [123,50,148]./256];
for odor_n = 1:3
    curr_color = color_vecs(odor_n, :);
    curr_trs = find(curr_tr_list == (odor_n - 1));
    curr_traces = curr_data(:, curr_trs, :);

    curr_means = squeeze(mean(curr_traces, 2, 'omitnan'));
    curr_ses = squeeze(std(curr_traces, [], 2, 'omitnan'))./sqrt(size(curr_traces, 2));

    figure(1)
    shadedErrorBar([], curr_means(:, 1), curr_ses(:, 1), {'Color', curr_color}, 1);
    hold on

    figure(2)
    shadedErrorBar([], curr_means(:, 2), curr_ses(:, 2), {'Color', curr_color}, 1);
    hold on

    figure(3)
    shadedErrorBar([], curr_means(:, 3), curr_ses(:, 3), {'Color', curr_color}, 1);
    hold on

end

bar_color = [0.8, 0.5, 0.5];

figure(1)
ylabel('probability(PA)');
ax_vals = axis;
plot([0, ax_vals(2)], [0.5, 0.5], '--', 'Color', bar_color);
set_xlabels_time(1, frame_time, 10);
fig_wrapup(1, [], [25, 30]);
add_stim_bar(1, stim_frs, bar_color);

figure(2)
ylabel('probability(BA)');
ax_vals = axis;
plot([0, ax_vals(2)], [0.5, 0.5], '--', 'Color', bar_color);
set_xlabels_time(1, frame_time, 10);
fig_wrapup(2, [], [25, 30]);
add_stim_bar(2, stim_frs, bar_color);


figure(3)
ylabel('probability(EL)');
ax_vals = axis;
plot([0, ax_vals(2)], [0.5, 0.5], '--', 'Color', bar_color);
set_xlabels_time(1, frame_time, 10);
fig_wrapup(3, [], [25, 30]);
add_stim_bar(3, stim_frs, bar_color);

end




