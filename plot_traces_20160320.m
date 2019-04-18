function [h, max_val] = plot_traces_20160320(dff_data_mat, t_list, stim_frame, stim_end_frame, stim_color, cell_n, odor_ni, fig_n, subplot_vec, normf, title_st)
%syntax: [h, max_val] = plot_traces_20160320(dff_data_mat, t_list, stim_frame, stim_end_frame, stim_color, cell_n, odor_ni, fig_n, subplot_vec, normf, title_st)
%This function plots single trial dF/F traces in grey with the averaged
%trace in a thicker black. It plots traces for a single cell, for a chosen
%odor number or across all odors.
%Mehrab Modi, 10/3/2014

if nargin == 4
    odor_n = 0;
    fig_n = 1;
    subplot_vec = [];
    normf = 0;
elseif nargin == 5
    fig_n = 1;
    subplot_vec = [];
    normf = 0;
elseif nargin == 6
    subplot_vec = [];
    normf = 0;
elseif nargin == 7
    normf = 0;
end


%close(findobj('type','figure','name','dF/F Traces'))


cell_traces = squeeze(dff_data_mat(:, cell_n, t_list));
ave_trace = nanmean(cell_traces, 2);
% cell_traces = cell_traces( 1:(size(dff_data_mat, 1)), : );
%ave_trace = ave_trace((stim_frame - 10):(size(dff_data_mat, 1)));

max_val = max(ave_trace);

if normf == 1
    cell_traces = cell_traces./max(ave_trace);
    ave_trace = ave_trace./max(ave_trace);
    
else
end


if isempty(subplot_vec) == 1
    figure('name', ['dF/F Traces' int2str(fig_n)]);
    grey = [0.6, 0.6, 0.6];
    h = plot(cell_traces);
    set(h, 'Color', grey);
    hold on
    plot(ave_trace, 'k', 'LineWidth', 3);
else
    figure(fig_n)
    subplot(subplot_vec(1), subplot_vec(2), subplot_vec(3));
    grey = [0.75, 0.75, 0.75];
    h = plot(cell_traces);
    set(h, 'Color', grey);
    hold on
    plot(ave_trace, 'k', 'LineWidth', 1);

end
%adding stimulus period patch
hold on
ax_lims = axis;
y_vec = [ax_lims(3:4)];
x_vec = [stim_frame, stim_end_frame];
y_vec = [y_vec, y_vec(2), y_vec(1)];
x_vec = [x_vec(1), x_vec(1), x_vec(2), x_vec(2)];
p = patch(x_vec', y_vec', stim_color);
set(p, 'FaceAlpha', 0.25);
set(p, 'EdgeColor', 'none');
set(gcf, 'Color', 'w')

title(title_st);
%keyboard
end