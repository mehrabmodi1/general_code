function fig_h = scattered_dot_plot(mat, fig_n, col_width, col_spacing, markersize, markercolor, marker_filled, with_lines, linecolor, xlabels, plot_mean, mean_color)
%syntax: fig_h = scattered_dot_plot(mat, fig_n, col_width, col_spacing, markersize, markercolor, marker_filled, with_lines, linecolor, xlabels, plot_mean, mean_color)
%This function plots the values in each column of mat as dots separated
%with a random scatter of width col_width and inter-column spacing as
%specified. Line spec can be used to specify marker style, color and size.
%If with_lines = 1, each point in a row will be plotted with the same
%random offset from its col center, with a line connecting these points
%across cols.

% col_width = 1;
% col_spacing = 1;
% mat = rand(100, 3);
% fig_n = 1;
% markersize = 5;
% markercolor = [0.5, 0.5, 0.5];

n_cols = size(mat, 2);

%allowing for explicitly defined colors for each column, or just one color
%for all
if size(markercolor, 1) == 1
    markercolor = repmat(markercolor, n_cols, 1);
else
end

fig_h = figure(fig_n);
saved_col_centers = zeros(1, n_cols);
if with_lines == 0
    for col_n = 1:n_cols
        curr_vec = mat(:, col_n);
        col_center = ( (col_n-1).*(col_width + (col_spacing)) + 0.5) + col_width./2;
        saved_col_centers(col_n) = col_center;
        r_vec = rand(length(curr_vec), 1).*col_width + ( (col_n-1).*(col_width + (col_spacing)) + 0.5);

        if marker_filled == 0
            plot(r_vec, curr_vec, 'O', 'markerSize', markersize, 'markerEdgeColor', markercolor(col_n, :))
        elseif marker_filled == 1
            plot(r_vec, curr_vec, 'O', 'markerSize', markersize, 'markerEdgeColor', markercolor(col_n, :), 'markerFaceColor', (markercolor(col_n, :)) )
        else
        end
        hold on
    end

    
elseif with_lines == 1
    %generating random offsets
    r_vec = rand(size(mat, 1), 1);
    for col_n = 1:n_cols
        col_center = ( (col_n-1).*(col_width + (col_spacing)) + 0.5) + col_width./2;
        saved_col_centers(col_n) = col_center;
        r_vecs(:, col_n) = r_vec.*col_width + ( (col_n-1).*(col_width + (col_spacing)) + 0.5);
    end
    
    %plotting each row as a line
    for row_n = 1:size(mat, 1)
        curr_row = mat(row_n, :);
        plot(r_vecs(row_n, :), curr_row, '-O', 'markerSize', markersize, 'markerEdgeColor', markercolor, 'Color', linecolor, 'lineWidth', 1)
        hold on
    end
    for row_n = 1:size(mat, 1)
        curr_row = mat(row_n, :);
        plot(r_vecs(row_n, :), curr_row, 'O', 'markerSize', markersize, 'markerEdgeColor', markercolor)
        hold on
    end
    
end

%plotting mean marker
if plot_mean == 1
   for col_n = 1:n_cols
        col_center = ( (col_n-1).*(col_width + (col_spacing)) + 0.5) + col_width./2;
        saved_col_centers(col_n) = col_center;
        curr_mean = mean(mat(:, col_n), 1, 'omitnan');
        curr_se = std(mat(:, col_n), 0, 1, 'omitnan')./sqrt(size(mat, 1) - sum(isnan(mat(:, col_n))));
        errorbar(col_center, curr_mean, curr_se, 'O', 'markerSize', markersize, 'markerEdgeColor', mean_color, 'markerFaceColor', mean_color, 'Color', mean_color, 'lineWidth', 2)
        
   end
    
  
    
else
end


ax = gca;
ax.XTick = saved_col_centers;
ax.XTickLabels = xlabels;
hold off
end