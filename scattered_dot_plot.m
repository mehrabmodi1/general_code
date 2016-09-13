function fig_h = scattered_dot_plot(mat, fig_n, col_width, col_spacing, markersize, markercolor, with_lines, linecolor, xlabels)
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

fig_h = figure(fig_n);
saved_col_centers = zeros(1, n_cols);
if with_lines == 0
    for col_n = 1:n_cols
        curr_vec = mat(:, col_n);
        col_center = ( (col_n-1).*(col_width + (col_spacing)) + 0.5) + col_width./2;
        saved_col_centers(col_n) = col_center;
        r_vec = rand(length(curr_vec), 1).*col_width + ( (col_n-1).*(col_width + (col_spacing)) + 0.5);


        plot(r_vec, curr_vec, 'O', 'markerSize', markersize, 'markerEdgeColor', markercolor)

        hold on
    end

    hold off
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
    hold off
end

ax = gca;
ax.XTick = saved_col_centers;
ax.XTickLabels = xlabels;

end