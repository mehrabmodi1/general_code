function [] = fig_wrapup(fig_n, script_name, varargin)
%Syntax: [] = fig_wrapup(fig_n, script_name)
%This function changes figure window size, gets rid of top and right side
%grid lines, adjusts ticklabel font size, and changes size, orientation of
%axis ticks. script_name is a string containing the name of the script
%calling fig_wrapup (obtained with mfilename).

% %testing lines
% vec = rand(1, 100);
% fig_h = figure(1);
% fig_n = 1;
% plot(vec)

plot_height = 100;
plot_width = 125;
%allowing custom-specified plot height and width
if isempty(varargin) == 0
    plot_height = varargin{1}(1);
    plot_width = varargin{1}(2);
else
end

axis_font_size = 7.5;
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
ax_handle.LineWidth = 1.2;
set(gca, 'fontname', 'Arial');

%adding script name to figure
if isempty(script_name) == 0
    add_scriptname(fig_n, script_name);
else
end

end