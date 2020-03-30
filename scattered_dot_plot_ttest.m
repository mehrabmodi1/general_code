function fig_h = scattered_dot_plot_ttest(mat, fig_n, col_width, col_spacing, markersize, markercolor, marker_filled, with_lines, linecolor, xlabels, plot_mean, mean_color, st_test_type)
%syntax: fig_h = scattered_dot_plot_ttest(mat, fig_n, col_width, col_spacing, markersize, markercolor, marker_filled, with_lines, linecolor, xlabels, plot_mean, mean_color, st_test_type)
%This function plots the values in each column of mat as dots separated
%with a random scatter of width col_width and inter-column spacing as
%specified. Line spec can be used to specify marker style, color and size.
%With_lines if not [], is a list of pairs of columns to be connected by
%lines.

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

if st_test_type == 1        %case with paired sample comparisons
   [p_vec_out] = stat_testing(mat, st_test_type); 
elseif st_test_tupe == 2    %case with two sample comparisons
    [p_vec_out] = stat_testing(mat, st_test_type);
else
end

fig_h = figure(fig_n);
saved_col_centers = zeros(1, n_cols);
if isempty(with_lines) == 1
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
        
        %adding p val label
        if rem(col_n, 2) == 0
            curr_veci = [mat(:, col_n); mat(:, (col_n - 1))];
            curr_p_val = round(p_vec_out(col_n./2, 1), 3);
            if curr_p_val ~= 0
                p_label = ['p = ', num2str(curr_p_val)];
            elseif curr_p_val == 0
                p_label = ['p < ', '0.001'];
            else
            end
            x_pos = mean(r_vec) - col_spacing;
            y_pos = max(curr_vec) + (max(max(mat)).*0.2);
            text(x_pos, y_pos, p_label);
        else
        end
    end

    
elseif isempty(with_lines) == 0
    %generating random offsets
    r_vec = rand(size(mat, 1), 1);
    for col_n = 1:n_cols
        col_center = ( (col_n-1).*(col_width + (col_spacing)) + 0.5) + col_width./2;
        saved_col_centers(col_n) = col_center;
        r_vecs(:, col_n) = r_vec.*col_width + ( (col_n-1).*(col_width + (col_spacing)) + 0.5);
                
        %adding p_val label
        if rem(col_n, 2) == 0
            curr_vec = [mat(:, col_n); mat(:, (col_n - 1))];
            curr_p_val = round(p_vec_out(col_n./2, 1), 3);
            if curr_p_val ~= 0
                p_label = ['p = ', num2str(curr_p_val)];
            elseif curr_p_val == 0
                p_label = ['p < ', '0.001'];
            else
            end
            x_pos = mean(r_vecs(:, col_n)) - col_spacing;
            y_pos = max(curr_vec) + (max(max(mat)).*0.2);
            text(x_pos, y_pos, p_label);
            hold on
        else
        end
    end
    
    %plotting lines for pairs of cols as specified in with_lines
    for row_n = 1:size(mat, 1)
        for conn_pair_n = 1:size(with_lines, 1)
            curr_pair = with_lines(conn_pair_n, :);
            try
                curr_row = mat(row_n, curr_pair(1):curr_pair(2) );
           
                
            if marker_filled == 0
                plot(r_vecs(row_n, curr_pair(1):curr_pair(2)), curr_row, '-O', 'markerSize', markersize, 'markerEdgeColor', markercolor(curr_pair(1), :), 'Color', linecolor(curr_pair(1), :), 'lineWidth', 1)
                hold on
                plot(r_vecs(row_n, curr_pair(2)), curr_row(2), '-O', 'markerSize', markersize, 'markerEdgeColor', markercolor(curr_pair(2), :), 'Color', linecolor(curr_pair(1), :), 'lineWidth', 1)
            elseif marker_filled == 1
                plot(r_vecs(row_n, curr_pair(1):curr_pair(2)), curr_row, '-O', 'markerSize', markersize, 'markerEdgeColor', markercolor(curr_pair(1), :), 'markerFaceColor', markercolor(curr_pair(1), :), 'Color', linecolor(curr_pair(1), :), 'lineWidth', 1)
                hold on
                plot(r_vecs(row_n, curr_pair(2)), curr_row(2), '-O', 'markerSize', markersize, 'markerEdgeColor', markercolor(curr_pair(2), :), 'markerFaceColor', markercolor(curr_pair(2), :), 'Color', linecolor(curr_pair(1), :), 'lineWidth', 1)
                
            else
            end
            
            catch
               keyboard
            end
            hold on
        end
    end
%     for row_n = 1:size(mat, 1)
%         curr_row = mat(row_n, :);
%         plot(r_vecs(row_n, :), curr_row, 'O', 'markerSize', markersize, 'markerEdgeColor', markercolor(col_n, :))
%         hold on
%     end
    
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
function [p_vec_out] = stat_testing(plot_mat, test_type)
    n_col_pairs = size(plot_mat, 2)./2;
    p_vec_out = zeros(n_col_pairs, 1) + nan;
    for col_pair_n = 0:(n_col_pairs - 1)
        col1 = (col_pair_n.*2) + 1;
        col2 = col1 + 1;
        if test_type == 1
            [del, p] = ttest(plot_mat(:, col1), plot_mat(:, col2));
        elseif test_type == 2
            [del, p] = ttest2(plot_mat(:, col1), plot_mat(:, col2));
        else
        end
        try
            p_vec_out((col_pair_n + 1), 1) = p;
        catch
            keyboard
    end
    

end

end