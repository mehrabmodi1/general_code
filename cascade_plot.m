function cascade_plot(figure_n, data_mat, linecolour, linewidth, x_offset_multiplier, y_offset_multiplier, ypatch_height, max_y_ax)
%Syntax: cascade_plot(figure_n, data_mat, linecolour, linewidth, x_offset_multiplier, y_offset_multiplier, ypatch_height, max_y_ax)
%This function creates a cascade plot of multiple traces, offset in y, and
%x to make all the traces visible. Line colour can be specified. 
%data_mat dims should be time, cells, trials.

figure(figure_n)

%computing max y range of all the vectors
max_vals = max(data_mat, [], 1, 'omitnan');
min_vals = min(data_mat, [], 1, 'omitnan');
yrange = max(max_vals - min_vals);
sing_vec = mean(data_mat, 2, 'omitnan');
n_nans = sum(isnan(sing_vec));
x_range = size(data_mat, 1) - n_nans;

if isempty(max_y_ax) ~= 1
    yrange = max_y_ax;
else
end
y_offset = yrange.*y_offset_multiplier;

if sign(x_offset_multiplier) ~= 0
    x_offset = round(x_range.*x_offset_multiplier);
    x_pad = zeros((abs(x_offset).*size(data_mat, 2)), size(data_mat, 2)) + nan;

    if sign(x_offset_multiplier) == 1
        data_mat = [data_mat; x_pad];
    elseif sign(x_offset_multiplier) == -1
        data_mat = [x_pad; data_mat];
    else
    end
else
end

for vec_n = 1:size(data_mat, 2)
    if sign(x_offset_multiplier) ~= 0
        data_mat(:, vec_n) = circshift(data_mat(:, vec_n), ((vec_n - 1).*x_offset ), 1);
    else
    end
    data_mat(:, vec_n) = data_mat(:, vec_n) + ((vec_n - 1) .* y_offset);
end

plot(data_mat, 'Color', linecolour, 'lineWidth', linewidth);


grid off
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'ycolor', 'none')

ax_vals = axis;
% ax_vals(5) = ax_vals(5) + (ax_vals(6) - ax_vals(5) )./5; %getting rid of lines at the bottom of patches
% ax_vals(6) = max_y_ax;
% ax_vals(1) = ax_vals(1) + 2;
% ax_vals(2) = data_vec_length - 2;
axis(ax_vals);

%adding y-scale patch
a = gca;
orig_ax = axis;
new_ax = [4.8, 15, orig_ax(3:4)];
orig_pos = a.Position;
new_pos = [(orig_pos(1) - (orig_pos(2).*0.25) ), orig_pos(1), orig_pos(3), orig_pos(4)];
ax2 = axes('Position', new_pos);    %created new axes object in same figure
axis(new_ax);
set(ax2, 'Color', 'none', 'Visible', 'off');

x1 = .1;
x2 = 5;

%patch width
y1 = ax_vals(3);
y2 = ypatch_height + ax_vals(3);

%drawing patch
y_vec = [y1, y2, y2, y1];
x_vec = [x1, x1, x2, x2];

p = patch(x_vec', y_vec', [0, 0, 0] );

p.EdgeColor = [1, 1, 1];

axes(a)