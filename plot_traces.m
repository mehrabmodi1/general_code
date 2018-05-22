function [h, max_val] = plot_traces(dff_data_mat, odor_t_list, stim_frame, cell_n, odor_n, fig_n, subplot_vec, norm)
%syntax: [h, max_val] = plot_traces(dff_data_mat, odor_t_list, stim_frame, cell_n, odor_n, fig_n, subplot_vec, norm)
%This function plots single trial dF/F traces in grey with the averaged
%trace in a thicker black. It plots traces for a single cell, for a chosen
%odor number or across all odors.
%Mehrab Modi, 10/3/2014

% if nargin == 4
%     odor_n = 0;
%     fig_n = 1;
%     subplot_vec = [];
%     norm = 0;
% elseif nargin == 5
%     fig_n = 1;
%     subplot_vec = [];
%     norm = 0;
% elseif nargin == 6
%     subplot_vec = [];
%     norm = 0;
% elseif nargin == 7
%     norm = 0;
% end


close(findobj('type','figure','name','dF/F Traces'))

%building list of trials for chosen odour
if odor_n == 0
    t_list = 1:size(dff_data_mat, 3);
else
end
   
if odor_n > 0
    t_list = find(odor_t_list == odor_n);
else
end


cell_traces = squeeze(dff_data_mat(:, cell_n, t_list));
ave_trace = nanmean(cell_traces, 2);

cell_traces = cell_traces( (stim_frame - 10):(size(dff_data_mat, 1)), : );
ave_trace = ave_trace((stim_frame - 10):(size(dff_data_mat, 1)));

max_val = max(ave_trace);

if norm == 1
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
    plot(ave_trace, 'k', 'LineWidth', 3);

end
set(gcf, 'Color', 'w')

end