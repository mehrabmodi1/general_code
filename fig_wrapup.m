function [] = fig_wrapup(fig_n)
%Syntax: [] = fig_wrapup(fig_n)
%This function changes figure window size, gets rid of top and right side
%grid lines, adjusts ticklabel font size, and changes size, orientation of
%axis ticks.

% %testing lines
% vec = rand(1, 100);
% fig_h = figure(1);
% fig_n = 1;
% plot(vec)


plot_height = 200;
plot_width = 280;
axis_font_size = 15;
tick_length = [0.02, 0.05];
line_widths = 2;
marker_sizes = 1;
marker_sizes_f = 4;

fig_h = figure(fig_n);
set(fig_h, 'Position', [100, 100, 100 + plot_width, 100 + plot_height]);
ax_handle = gca;
ax_handle.FontSize = axis_font_size;
ax_handle.TickLength = [0.005, 0.005];
ax_handle.Box = 'off';
ax_handle.TickDir = 'out';

end