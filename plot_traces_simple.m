function [h, max_val] = plot_traces_simple(fig_n, dff_data_mat, cell_n, odor_n, dur_n, stim_mat_simple, stim_mat, frame_time, color_vec, norm, plot_mean)
%syntax: function [h, max_val] = plot_traces_simple(fig_n, dff_data_mat, cell_n, odor_n, dur_n, stim_mat_simple, stim_mat, frame_time, norm, plot_mean)
%This function plots single trial dF/F traces in grey with the averaged
%trace in a thicker black. It plots traces for a single cell, for a chosen
%odor number or across all odors.
%Mehrab Modi, 20180522

odor_list = unique(stim_mat_simple(:, 2));
odor_dur_list = unique(stim_mat_simple(:, 3));
curr_dur = odor_dur_list(dur_n);
odor_ni = odor_list(odor_n);

curr_trs = find(stim_mat_simple(:, 2) == odor_ni & stim_mat_simple(:, 3) == curr_dur);
stim_frs = compute_stim_frs(stim_mat, curr_trs(1), frame_time);

cell_traces = squeeze(dff_data_mat(:, cell_n, curr_trs));
ave_trace = mean(cell_traces, 2, 'omitnan');

max_val = max(ave_trace);

if norm == 1
    cell_traces = cell_traces./max(ave_trace);
    ave_trace = ave_trace./max(ave_trace);
    
else
end

grey = [0.6, 0.6, 0.6];
figure(fig_n)
h = plot(cell_traces);
set(h, 'Color', grey);
if plot_mean == 1
    hold on
    plot(ave_trace, 'k', 'LineWidth', 3);
    hold off
else
end
set(gcf, 'Color', 'w')
set_xlabels_time(fig_n, frame_time, 5);
add_stim_bar(fig_n, stim_frs, color_vec(odor_n, :))

end