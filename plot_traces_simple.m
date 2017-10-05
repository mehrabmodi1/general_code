function [h, max_val] = plot_traces_simple(trace_mat, time_vec, norm)
%This function plots single trial dF/F traces in grey with the averaged
%trace in a thicker black. It plots traces for a single cell, for a chosen
%odor number or across all odors.
%Mehrab Modi, 10/3/2014

close(findobj('type','figure','name','dF/F Traces'))

cell_traces = trace_mat;
ave_trace = mean(cell_traces, 2, 'omitnan');

max_val = max(ave_trace);

if norm == 1
    cell_traces = cell_traces./max(ave_trace);
    ave_trace = ave_trace./max(ave_trace);
    
else
end

grey = [0.6, 0.6, 0.6];
h = plot(time_vec, cell_traces);
set(h, 'Color', grey);
hold on
plot(time_vec, ave_trace, 'k', 'LineWidth', 3);

set(gcf, 'Color', 'w')

end