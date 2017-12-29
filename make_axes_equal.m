function [] = make_axes_equal(fig_n, equality_line)
%Syntax: [] = make_axes_equal(fig_n, equality_line)
%This function sets the two axes of a plot to be equal. If equality_line is 1,
%it also plots a dotted, grey, equality line.

%test variables
% fig_n = 1;
% figure(1)
% plot(rand(1, 5))
% equality_line = 1;

figure(fig_n)
ax_vals = axis;

max_val = max([ax_vals(1, 2), ax_vals(1, 4)]);

ax_vals(1, 2) = max_val;
ax_vals(1, 4) = max_val;

%plotting equality line if specified
if equality_line == 1
    hold on
    plot([0, max_val], [0, max_val], '--', 'Color', [0.65, 0.65, 0.65], 'LineWidth', 2)
    hold off
else
end


end