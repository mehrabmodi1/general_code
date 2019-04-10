function [] = fig_wrapup_nonums(fig_n, script_name, dims)
%Syntax: [] = fig_wrapup_nonums(fig_n, script_name)
%This function changes figure window size, gets rid of top and right side
%grid lines, adjusts ticklabel font size, and changes size, orientation of
%axis ticks. script_name is a string containing the name of the script
%calling fig_wrapup (obtained with mfilename).

% %testing lines
% vec = rand(1, 100);
% fig_h = figure(1);
% fig_n = 1;
% plot(vec)


plot_height = 200;
plot_width = 250;
axis_font_size = 15;
tick_length = [0.02, 0.05];
line_widths = 2;
marker_sizes = 1;
marker_sizes_f = 4;

fig_h = figure(fig_n);
set(fig_h, 'Position', [100, 100, 100 + plot_width, 100 + plot_height]);
ax_handle = gca;
if isempty(intersect(dims, 1)) == 0
    ax_handle.XTickLabel = [];
elseif isempty(intersect(dims, 2)) == 0
    ax_handle.YTickLabel = [];
else
end

ax_handle.FontSize = axis_font_size;
ax_handle.TickLength = [0.005, 0.005];
ax_handle.Box = 'off';
ax_handle.TickDir = 'out';
ax_handle.LineWidth = 1.2;
%adding script name to figure
add_scriptname(fig_n, script_name);

end