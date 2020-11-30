function cascade_plot(figure_n, data_mat, linecolour, linewidth, x_offset_multiplier, y_offset_multiplier, ypatch_height, max_y_ax, highlight_pks, indicate_zeros)
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

%highlighting peaks with red circles
if highlight_pks == 1
    [pk_vals, pk_locs] = max(data_mat, [], 1, 'omitnan');
    hold on
    plot(pk_locs, pk_vals, 'rO');
else
end

%indicating 0 for each trace with a dotted line
if indicate_zeros == 1
    zeros_vec = (0:(size(data_mat, 2)) - 1) .* y_offset;
    zeros_mat = repmat(zeros_vec, size(data_mat, 1), 1);
    hold on
    plot(zeros_mat, ':', 'Color', [0.6, 0.6, 0.6], 'lineWidth', 0.75)
    hold off
else
end


grid off
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'ycolor', 'none')

ax_vals = axis;

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