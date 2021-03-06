function [fig_h] = bar_line_plot(fig_n, data_mat, line_color, pre_color, post_color, bar_width, bar_space)
%This function plots bars with error bars, overlaid with lines representing
%each paired sample that makes up the two bars.

fig_h = figure(fig_n);
n_pairs = size(data_mat, 2)./2;
mean_vecs = mean(data_mat, 1);
mean_vecs = reshape(mean_vecs, 2, []);
se_vecs = std(data_mat, [], 1, 'omitnan')./sqrt(size(data_mat, 1));
se_vecs = reshape(se_vecs, 2, []);


%plotting bars
b = bar(mean_vecs', 0.7, 'FaceColor', 'flat');
b(1).CData(1:n_pairs, :) = repmat(pre_color, n_pairs, 1);
b(2).CData(1:n_pairs, :) = repmat(post_color, n_pairs, 1);

hold on


%plotting errorbars 
for bar_n = 1:2
    xData = b(bar_n).XData + b(bar_n).XOffset;
    errorbar(xData, mean_vecs(bar_n, :), se_vecs(bar_n, :), 'Linestyle', 'none', 'CapSize', 8, 'Color', 'k', 'LineWidth', 1);
    xData_mat(bar_n, :) = xData;
end

%plotting paired obs lines
for pair_n = 1:n_pairs
    x_vals = xData_mat(:, pair_n)';
    x_vals = repmat(x_vals, size(data_mat, 1), 1 );
    data_col_n = (pair_n - 1).*2 + 1;
    y_vals = data_mat(:, data_col_n:(data_col_n + 1));
    
    plot(x_vals', y_vals', 'Color', line_color, 'LineWidth', 1)
end

hold off
