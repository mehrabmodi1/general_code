function [resp_areas_sorted, fig_h, c_vals_mat] = data_quality_resp_matrix(resp_areas, stim_mat_simple, sig_cell_mat, fig_n, plot_matrix)
%Syntax: [resp_areas_sorted, fig_h] = data_quality_resp_matrix(resp_area_mat, stim_mat_simple, plot_matrix)
%This function re-arranges response areas trial by trial, grouping trials
%with the same odor stimulus. This lets the user check data quality at a glance.


odor_list = unique(stim_mat_simple(:, 2));
odor_dur_list = unique(stim_mat_simple(: ,3));

resp_areas_sorted = [];
for odor_n = 1:length(odor_list)
    odor_ni = odor_list(odor_n);
    for odor_dur_n = 1:length(odor_dur_list)
        curr_dur = odor_dur_list(odor_dur_n);
        curr_trs = find(stim_mat_simple(:, 2) == odor_ni & stim_mat_simple(:, 3) == curr_dur);
        resp_areas_sorted = [resp_areas_sorted, resp_areas(:, curr_trs)]; 
        mean_resp_areas(:, odor_n) = mean(resp_areas(:, curr_trs), 2, 'omitnan');
        curr_sig_cells = sig_cell_mat(odor_ni, odor_dur_n);
        curr_insig_cells = find(curr_sig_cells == 0);
        mean_resp_areas(curr_insig_cells) = 0;
        del = find(mean_resp_areas > 2);
        mean_resp_areas(del) = 2;
    end
end

c = corrcoef(mean_resp_areas, 'Rows', 'complete');
c_vals_mat = [c(1, 2), c(1, 3), c(2, 3)];

if plot_matrix == 1
    fig_h = figure(fig_n);
    imagesc(resp_areas_sorted', [0, 5])
    xlabel('cell number')
    ylabel('trials')
else
    fig_h = [];
end