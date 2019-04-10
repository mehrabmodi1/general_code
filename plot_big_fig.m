function [] = plot_big_fig(fig_n)
%Syntax: [] = plot_big_fig(fig_n)
%This function identifies the biggest monitor available and resizes fig_n
%to be as big as possible on this monitor

% figure(1)
% plot(rand(5, 5))
% fig_n = 1;


fig_size = get(gcf, 'Position');
fig_size = fig_size(3:4);

figure(fig_n)
monitor_pos = get(0,'MonitorPositions');
[small_dim, big_monitor_n] = max(min(monitor_pos(:, 3:4), [], 2));
mon_size = monitor_pos(big_monitor_n, 3:4);

conversion_ratio = min(mon_size./fig_size);     %ratio of size of figure to size of monitor on each dimension
fig_size_n = floor(fig_size.*conversion_ratio.*0.85); %re-calculating size of figure, multiplying by smaller monitor/fig ratio

origin = monitor_pos(big_monitor_n, 1:2);
origin_n = [(origin(1)+ (mon_size(1)./2 - fig_size_n(1)./2)), origin(2)];

set(gcf, 'Position', [origin_n, fig_size_n]);

