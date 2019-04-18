function [] = cascade_plot(fig_n, data_mat, plot_vs, ydim, norm, offset_mult, linew, line_colour)
%cascade_plot syntax: cascade_plot(data_mat, ydim, offset, linew) cascade_plot plots 
%multiple traces in one window with a y-displacement for each row. 
%Each row is individually normalised to its maximum value. cascade_plot can 
%only recieve a 2-D matrix as input, and will consider 
%ydim as the y axis for the cascade plot. 'offset' is a multiplicative
%factor that changes the y-offset between traces. 'linew' sets the
%thickness of the lines used to make the plot.

%testing inputs
% fig_n = 1;
% data_mat = rand(100, 10);
% plot_vs = 1:1:100;
% ydim = 2;
% norm = 0;
% offset_mult = 1.1;
% linew = 2;
% line_colour = [.6, .6, .6];


if length(size(data_mat)) > 2
    error('Data_mat has > 2 dimensions. Only 2-D matrix input to cascade_plot.')
else
end

curr_handle = figure(fig_n);

if ydim == 2
    data_mat = data_mat';
else
end

max_max_val = max(max(data_mat));

for row_no = 1:size(data_mat, 1)
    max_val = max(data_mat(row_no, :));
    if norm == 1
        data_mat(row_no, :) = data_mat(row_no, :)./max_val + row_no.*1.2.*offset_mult;
    elseif norm ~= 1
        data_mat(row_no, :) = data_mat(row_no, :) + row_no.*max_max_val.*offset_mult;
    end
end

plot(plot_vs, data_mat', 'LineWidth', linew, 'Color', line_colour)

%wrapping up figure and getting rid of the y axis markers
disp('Dont call fig_wrapup after this function.')
plot_height = 200;
plot_width = 280;
axis_font_size = 15;
line_widths = 2;
marker_sizes = 1;
marker_sizes_f = 4;

fig_h = figure(fig_n);
set(fig_h, 'Position', [100, 100, 100 + plot_width, 100 + plot_height]);
ax_handle = gca;
ax_handle.FontSize = axis_font_size;
ax_handle.TickLength = [0.005, 0];
ax_handle.Box = 'off';
ax_handle.TickDir = 'out';
ytickmin = get(ax_handle, 'YTick');
ytickmin = ytickmin(2);
disp(['y-bar height: ' int2str(ytickmin)])
set(gca,'YTick',[], 'YColor', [1, 1, 1])
    
%adding patch of height ytickmin
figure(fig_n);
a = gca;
orig_ax = axis;
new_ax = [0, 0.3, 0, orig_ax(4)];
orig_pos = a.Position;
new_pos = [(orig_pos(1) - orig_pos(1).*0.5), orig_pos(2), .05, orig_pos(4) ];
ax2 = axes('Position', new_pos);    %created new axes object in same figure
axis(new_ax);
set(ax2, 'Color', 'none', 'Visible', 'off');

%patch height
y1 = 0;
y2 = ytickmin;
disp(['bar height = ' num2str(ytickmin)])

%patch width
x1 = 0.1;
x2 = 0.2;

%drawing patch
y_vec = [y1, y2, y2, y1];
x_vec = [x1, x1, x2, x2];
p = patch(x_vec', y_vec', [0, 0, 0]);
p.EdgeColor = [1, 1, 1];

children = allchild(curr_handle);
axes(children(size(children, 1)));


                        
