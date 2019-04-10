clear all
close all

data = [-71.61, -56.74, -49.47, -43.215, -53.72, -47.13, -42.66, -60.285;...
        -9.785, -1.6, -21.53, -14.335, -22, -0.62, -0.925, -32.395];
data = data'./100;
col_width = .3;
col_space = 1;
fig1h = scattered_dot_plot(data, 1, col_space, col_width, 10, [.25, .25, .25], 0, 'k', {'60s odour', '1s odour'});
hold on
means = mean(data);
ses = std(data)./sqrt(size(data, 1));
x_vals = [col_space, (2.*col_space + col_width)];
errorbar(x_vals, means, ses, 'LineStyle', 'none', 'MarkerSize', 10, 'Color', 'r', 'LineWidth', 2)
plot(x_vals, means, '.r', 'MarkerSize', 30)
ax = axis;
plot([0, 100], [0, 0], '--r', 'LineWidth', 3)
axis(ax);
ylabel('avoidance index')
fig_wrapup(1);

[h, p] = ttest2(data(:, 1), data(:, 2))