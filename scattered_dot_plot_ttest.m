function [fig_h, r_vecs_saved, saved_col_centers] = scattered_dot_plot_ttest(mat, fig_n, col_width, col_spacing, markersize, markercolor, marker_filled, with_lines, linecolor, xlabels, plot_mean, mean_color, st_test_type, p_val, coh_d, varargin)
%syntax: fig_h = scattered_dot_plot_ttest(mat, fig_n, col_width, col_spacing, markersize, markercolor, marker_filled, with_lines, linecolor, xlabels, plot_mean, mean_color, st_test_type, p_val, coh_d, bar_thickness, force_means)
%This function plots the values in each column of mat as dots separated
%with a random scatter of width col_width and inter-column spacing as
%specified. Line spec can be used to specify marker style, color and size.
%With_lines if not [], is a list of pairs of columns to be connected by
%lines. st_test_type = 1, paired-samples. st_test_type = 2, independent
%samples.


%if using the beeswarm function to generate r_offsets, use these defaults
% mat = rand(100, 3);
% fig_n = 1;
% markersize = 5;
% markercolor = [0.5, 0.5, 0.5];

if isempty(varargin) == 0
    bar_thickness = varargin{1};
    force_means = varargin{2};
    if length(varargin) >= 3
        wrapup_vars = varargin{3};
        if isempty(wrapup_vars) == 1
            wrapup_vars{1} = [100, 120];
            wrapup_vars{2} = 0.6;
        else
        end
        mean_markersize = markersize;
    else
        wrapup_vars{1} = [100, 120];
        wrapup_vars{2} = 0.6;
        mean_markersize = markersize;
    end    
    
    if force_means == 'force_mean'
        force_means = 1;
    else
        force_means = 0;
    end
    
    if length(varargin)>= 4
        plot_pval = varargin{4};
    else
        plot_pval = 1;
        mean_markersize = markersize;
    end
    
    if length(varargin)>= 5
        mean_markersize = varargin{5};
    else
        mean_markersize = markersize;
    end
    
else
    bar_thickness = 2;
    force_means = 0;
    wrapup_vars{1} = [100, 120];
    wrapup_vars{2} = 0.6;
    plot_pval = 1;
    mean_markersize = markersize;
end


n_cols = size(mat, 2);
violin_cutoff = 500;    %if n points in any column crosses this thresh, all columns rendered as violins
if size(mat, 1) > violin_cutoff
    plot_violins = 1;
else
    plot_violins = 0;
end

%allowing for explicitly defined colors for each column, or just one color
%for all
if size(markercolor, 1) == 1
    markercolor = repmat(markercolor, n_cols, 1);
else
end


[p_vec_out, effect_sizes_out, min_difs_out] = stat_testing(mat, st_test_type, plot_mean, p_val); 

curr_class = whos('fig_n');
if strcmp(curr_class.class, 'double') == 1
    fig_h = figure(fig_n);
else
    fig_h = fig_n;
end

if isempty(with_lines) == 0
    bees_on = 0;
else
    bees_on = 1;
    col_width = 1;
    %col_spacing = 0;
end
   

if bees_on == 1
    %using beeswarm to generate random offset vectors that retain distribution shape for all data columns
    %setting up input vectors for beeswarm
    x_vec = zeros((size(mat, 1).*size(mat, 2)), 1);
    y_vec = x_vec;
    for col_n = 1:size(mat, 2)
        st_pt = (col_n - 1).*size(mat, 1) + 1;
        stp_pt = col_n.*size(mat, 1);
        x_vec(st_pt:stp_pt, 1) = col_n;             %generating vector of group numbers (col numbers)
        y_vec(st_pt:stp_pt, 1) = mat(:, col_n);     %generating a single vector of all the columns
    end
    r_vecs = beeswarm_noplot(x_vec, y_vec, wrapup_vars, 'sort_style','rand', 'corral_style', 'gutter');     %vector of x positions for each point
    r_mat = zeros(size(mat, 1), size(mat, 2));
    for col_n = 1:size(mat, 2)
        st_pt = (col_n - 1).*size(mat, 1) + 1;
        stp_pt = col_n.*size(mat, 1);
        r_mat(:, col_n) = r_vecs(st_pt:stp_pt, 1);  
    end
    
    %subtracting column offsets so each column ranges from 0 to 1
    offset_vec = [1:1:size(r_mat, 2)];
    r_mat = r_mat - repmat(offset_vec, size(r_mat, 1), 1);
    
    clear r_vecs
else
end

%allowing different col-spacings for within vs across pair columns
if length(col_spacing) == 1
    col_spacing_vec = [col_spacing, col_spacing];
else
    col_spacing_vec = col_spacing;
end

saved_col_centers = zeros(1, n_cols);
if isempty(with_lines) == 1
    x_vec = [];
    r_vecs_saved = [];
    
    for col_n = 1:n_cols
        if rem(col_n, 2) == 1
            col_spacing = col_spacing_vec(1);
        elseif rem(col_n, 2) == 0
            col_spacing = col_spacing_vec(2);
        else
        end
        curr_vec = mat(:, col_n);
        %col_center = ( (col_n-1).*(col_width + (col_spacing)) + 0.5) + col_width./2;
        if col_n > 1
            col_center = saved_col_centers(col_n-1) + col_spacing + (col_width./2);
        elseif col_n == 1
            col_center = 0 + col_spacing + (col_width./2);
        else
        end
        
        saved_col_centers(col_n) = col_center;
        
        if bees_on == 1
            r_vec = r_mat(:, col_n).*col_width + col_center;        %taking x positions for each point from beeswarm generated r_mat
        else
            %r_vec = rand(length(curr_vec), 1).*col_width + ( (col_n-1).*(col_width + (col_spacing)) + 0.5);
            r_vec = rand(length(curr_vec), 1).*col_width + col_center;
        end
        
        r_vecs_saved = [r_vecs_saved, r_vec];
        %r_vec_center = 0.5.*col_width + ( (col_n-1).*(col_width + (col_spacing)) + 0.5);    %for use with a violin plot, if needed
        r_vec_center = 0.5.*col_width + col_center;    %for use with a violin plot, if needed
        if plot_violins == 1
            x_vec = [x_vec, r_vec_center];
            hold on
        elseif plot_violins == 0
            if marker_filled == 0
                plot(r_vec, curr_vec, 'O', 'markerSize', markersize, 'markerEdgeColor', markercolor(col_n, :));
            elseif marker_filled == 1
                plot(r_vec, curr_vec, 'O', 'markerSize', markersize, 'markerEdgeColor', markercolor(col_n, :), 'markerFaceColor', (markercolor(col_n, :)) );
            else
            end
            hold on
        end
        
        %adding p val label
        if rem(col_n, 2) == 0
            curr_veci = [mat(:, col_n); mat(:, (col_n - 1))];
            curr_p_val = round(p_vec_out(col_n./2, 1), 3);
            curr_d = round(effect_sizes_out(col_n./2, 1), 3);
            if plot_mean == 1
                try
                    curr_dif = round(min_difs_out(col_n./2, 1), 3);
                catch
                    keyboard
                end
                curr_ac_dif = round(min_difs_out(col_n./2, 2), 3);
            elseif plot_mean == 2
                curr_dif = [];
                curr_ac_dif = [];
            else
            end
            if plot_mean == 1
                min_n = min_difs_out(col_n./2, 3);
                ac_n = min_difs_out(col_n./2, 4);
            elseif plot_mean == 2
                min_n = [];
                ac_n = [];
            else
            end
            
            if curr_p_val ~= 0
                p_label = ['p = ', num2str(curr_p_val)];
            elseif curr_p_val == 0
                p_label = ['p < ', '0.001'];
            else
            end
            if curr_d ~= 0
                d_label = ['Cohen''s d = ', num2str(curr_d)];
            elseif curr_d == 0
                d_label = ['Cohen''s d  < ', '0.001'];
            else
            end
            n_label = ['min n = ' num2str(min_n), ', n = ' num2str(ac_n)];
            if curr_dif ~= 0
                dif_label = ['min. dif = ', num2str(curr_dif), ', dif = ' num2str(curr_ac_dif)];
            elseif curr_dif == 0
                dif_label = ['min. dif  < ', '0.001, dif = ', num2str(curr_ac_dif)];
            elseif isempty(curr_dif) == 1
                dif_label = [];
                n_label = [];
            else
            end
            x_pos = mean(r_vec) - col_spacing.*1.2;
            y_pos = max(curr_vec) + (max(max(mat)).*0.1);
            
            if plot_pval == 1
                text(x_pos, y_pos, p_label, 'FontName', 'Arial', 'FontSize', 7.5);
                if coh_d == 1
                    text(x_pos, (y_pos + (max(max(mat)).*0.1)), d_label, 'FontName', 'Arial', 'FontSize', 7.5);
                else
                end
            else
            end
%             text(x_pos, (y_pos + 2.*(max(max(mat)).*0.1)), dif_label);
%             text(x_pos, (y_pos + 3.*(max(max(mat)).*0.1)), n_label);
        else
        end
    end
    
    if plot_violins == 1
        try

            violin(mat, 'x', x_vec, 'facecolor', markercolor, 'edgecolor','none', 'mc', [],'medc', [], 'width_multiplier', col_spacing, 'facealpha', 1);
            
        catch
            keyboard
        end
    else
    end
    
elseif isempty(with_lines) == 0 | with_lines == 0
    
    %generating random offsets
    r_vec = rand(size(mat, 1), 1);
    for col_n = 1:n_cols
        if rem(col_n, 2) == 1
            col_spacing = col_spacing_vec(1);
        elseif rem(col_n, 2) == 0
            col_spacing = col_spacing_vec(2);
        else
        end
        
        %col_center = ( (col_n-1).*(col_width + (col_spacing)) + 0.5) + col_width./2;
        if col_n > 1
            col_center = saved_col_centers(col_n - 1) + col_spacing + (col_width./2);
        elseif col_n == 1
            col_center = 0 + col_spacing + (col_width./2);
        else
        end
        
        saved_col_centers(col_n) = col_center;
        
        
        r_vecs(:, col_n) = r_vec.*col_width + col_center;
                
        %adding p_val label
        if rem(col_n, 2) == 0
            curr_vec = [mat(:, col_n); mat(:, (col_n - 1))];
            curr_p_val = round(p_vec_out(col_n./2, 1), 3);
            curr_d = round(effect_sizes_out(col_n./2, 1), 3);
            if plot_mean == 1
                curr_dif = round(min_difs_out(col_n./2, 1), 3);
                curr_ac_dif = round(min_difs_out(col_n./2, 2), 3);
                min_n = min_difs_out(col_n./2, 3);
                ac_n = min_difs_out(col_n./2, 4);
            else
                curr_dif = [];
                curr_ac_dif = [];
                min_n = [];
                ac_n = [];
            end
            
            if curr_p_val ~= 0
                p_label = ['p = ', num2str(curr_p_val)];
            elseif curr_p_val == 0
                p_label = ['p < ', '0.001'];
            else
            end
            if curr_d ~= 0
                d_label = ['Cohen''s d = ', num2str(curr_d)];
            elseif curr_d == 0
                d_label = ['Cohen''s d  < ', '0.001'];
            else
            end
            n_label = ['min n = ' num2str(min_n), ', n = ' num2str(ac_n)];
            if curr_dif ~= 0
                dif_label = ['min. dif = ', num2str(curr_dif), ', dif = ' num2str(curr_ac_dif)];
            elseif curr_dif == 0
                dif_label = ['min. dif  < ', '0.001, dif = num2str(curr_ac_dif)'];
            elseif isempty(curr_dif) == 1
                dif_label = [];
                n_label = [];
            else
            end
            x_pos = mean(r_vecs(:, col_n)) - col_spacing.*1.2;
            y_pos = max(curr_vec) + (max(max(mat)).*0.1);
            text(x_pos, y_pos, p_label, 'FontName', 'Arial', 'FontSize', 7.5);
            if coh_d == 1
                text(x_pos, (y_pos + (max(max(mat)).*0.1)), d_label, 'FontName', 'Arial', 'FontSize', 7.5);
            else
            end
%             text(x_pos, (y_pos + 2.*(max(max(mat)).*0.1)), dif_label);
%             text(x_pos, (y_pos + 3.*(max(max(mat)).*0.1)), n_label);
            hold on
            
        else
        end
    end
    r_vecs_saved = r_vecs;
    
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
if size(mean_color, 1) > 1
   mean_color_vec = mean_color;
else
   mean_color_vec = repmat(mean_color, n_cols, 1);
end
if plot_mean == 1
   for col_n = 1:n_cols
        mean_color = mean_color_vec(col_n, :);
        %col_center = ( (col_n-1).*(col_width + (col_spacing)) + 0.5) + col_width./2;
        col_center = saved_col_centers(col_n);
        %saved_col_centers(col_n) = col_center;
        curr_mean = mean(mat(:, col_n), 1, 'omitnan');
        curr_se = std(mat(:, col_n), 0, 1, 'omitnan')./sqrt(size(mat, 1) - sum(isnan(mat(:, col_n))));
        errorbar(col_center, curr_mean, curr_se, 'O', 'markerSize', mean_markersize, 'markerEdgeColor', mean_color, 'markerFaceColor', mean_color, 'Color', mean_color, 'lineWidth', bar_thickness)

   end
elseif plot_mean == 2        %plotting median and quantiles
    for col_n = 1:n_cols
        mean_color = mean_color_vec(col_n, :);
        %col_center = ( (col_n-1).*(col_width + (col_spacing)) + 0.5) + col_width./2;
        col_center = saved_col_centers(col_n);
        %saved_col_centers(col_n) = col_center;
        if force_means == 0
            curr_mean = median(mat(:, col_n), 1, 'omitnan');
            curr_se_up = quantile(mat(:, col_n), 0.75, 1) - curr_mean;      %subtracting mean because errorbar needs length of bar as input, not absolute y val
            curr_se_down = curr_mean - quantile(mat(:, col_n), 0.25, 1);    %%subtracting mean because errorbar needs length of bar as input, not absolute y val
            errorbar(col_center, curr_mean, curr_se_down, curr_se_up, 'O', 'markerSize', mean_markersize, 'markerEdgeColor', mean_color, 'markerFaceColor', mean_color, 'Color', mean_color, 'lineWidth', bar_thickness)
        elseif force_means == 1
            curr_mean = mean(mat(:, col_n), 1, 'omitnan');
            curr_se = std(mat(:, col_n), 0, 1, 'omitnan')./sqrt(size(mat, 1) - sum(isnan(mat(:, col_n))));
            errorbar(col_center, curr_mean, curr_se, 'O', 'markerSize', mean_markersize, 'markerEdgeColor', mean_color, 'markerFaceColor', mean_color, 'Color', mean_color, 'lineWidth', bar_thickness)
        else
        end
        
    end
else
    for col_n = 1:n_cols
        %col_center = ( (col_n-1).*(col_width + (col_spacing)) + 0.5) + col_width./2;
        col_center = saved_col_centers(col_n);
        %saved_col_centers(col_n) = col_center;
    end
end

hold on
ax_vals = axis();
plot([ax_vals(1), ax_vals(2)], [0, 0], ':', 'lineWidth', 2, 'Color', [0.6, 0.6, 0.6]);


ax = gca;
ax.XTick = saved_col_centers;
ax.XTickLabels = xlabels;
hold off


end


function [p_vec_out, effect_sizes_out, min_dif_out] = stat_testing(plot_mat, test_type, plot_mean, p_val)
    if rem(size(plot_mat, 2), 2) == 0
        
    else
        %odd number of columns, no pair-wise comparisons being made.
        p_vec_out = zeros(size(plot_mat, 2)) + nan;
        effect_sizes_out = p_vec_out;
        min_dif_out = zeros(size(plot_mat, 2), 4) + nan;
        return
    end
    n_col_pairs = size(plot_mat, 2)./2;
    
    p_vec_out = zeros(n_col_pairs, 1) + nan;
    effect_sizes_out = zeros(n_col_pairs, 1) + nan;
    for col_pair_n = 0:(n_col_pairs - 1)
        col1 = (col_pair_n.*2) + 1;
        col2 = col1 + 1;
        n = size(plot_mat, 1);
            
            if plot_mean == 1   %means, assume normal distributions
                if test_type == 1       %paired samples
                    [del, p] = ttest(plot_mat(:, col1), plot_mat(:, col2));
                    effect_size = computeCohen_d(plot_mat(:, col1), plot_mat(:, col2), 'paired');
                    paired_difs = plot_mat(:, col1) - plot_mat(:, col2);    %for paired-samples, sample distribution is distribution of paired differences.
                    mean1 = mean( paired_difs, 'omitnan');                  %actual mean of paired differences
                    sd1 = std(paired_difs, [], 'omitnan');
                    min_mean1 = sampsizepwr('t',[0 sd1],[],(1 - p_val.*2), n);        %null hypothesis mean is difference of 0
                    needed_n = sampsizepwr('t',[0 sd1], mean1, (1 - p_val.*2), []);
                    min_dif = abs(min_mean1 - 0); 
                    ac_dif = mean1;
                elseif test_type == 2   %independent samples
                    [del, p] = ttest2(plot_mat(:, col1), plot_mat(:, col2));
                    effect_size = computeCohen_d(plot_mat(:, col1), plot_mat(:, col2), 'independent');
                    mean1 = mean(plot_mat(:, col1), 'omitnan');
                    ac_mean2 = mean(plot_mat(:, col2), 'omitnan');
                    sd1 = std(plot_mat(:, col1), [], 'omitnan');
                    try
                        mean2 = sampsizepwr('t2',[mean1 sd1],[],(1 - p_val.*2), n);        %two-tailed test
                    catch
                        keyboard
                    end
                    needed_n = sampsizepwr('t',[0 sd1], mean1, (1 - p_val.*2), []);
                    min_dif = abs(mean2 - mean1);
                    ac_dif = abs(ac_mean2 - mean1);
                else
                end
                
            elseif plot_mean == 2   %medians, non-parametric tests
                if test_type == 1       %paired samples
                    p = signrank(plot_mat(:, col1), plot_mat(:, col2));
                    effect_size = computeCohen_d(plot_mat(:, col1), plot_mat(:, col2), 'paired');
                    min_dif = [];
                    ac_dif = [];
                    needed_n = [];
                    n = [];
                elseif test_type == 2   %independent samples
                    p = ranksum(plot_mat(:, col1), plot_mat(:, col2));    %Mann Whitney U test
                    effect_size = computeCohen_d(plot_mat(:, col1), plot_mat(:, col2), 'independent');
                    min_dif = [];
                    needed_n = [];
                    ac_dif = [];
                    n = [];
                else
                end
                
            else
            end
        try
            p_vec_out((col_pair_n + 1), 1) = p;
            effect_sizes_out((col_pair_n + 1), 1) = effect_size;
            if test_type == 1 && plot_mean == 1
                min_dif_out((col_pair_n + 1), 1:4) = [min_dif, ac_dif, needed_n, n];
            elseif test_type == 2 && plot_mean == 1
                min_dif_out((col_pair_n + 1), 1:4) = [min_dif, ac_dif, needed_n, n];
            else
                min_dif_out = [];
            end
        catch
            keyboard
        end
    

    end
    
    if test_type == 1   %paired-samples
        if plot_mean == 1   %use means, parametric test
            disp('Mean with SEs, paired-sample T test.')
        elseif plot_mean == 2  %use medians, non-parametric test
            disp('Median with quartiles, Wilcoxon signed rank test.')
        else
        end
    elseif test_type == 2   %independent samples
        if plot_mean == 1   %use means, parametric test
            disp('Mean with SEs, two-sample T test.')
        elseif plot_mean == 2  %use medians, non-parametric test
            disp('Median with quartiles, Mann Whitney U test.')
        else
        end
    else
    end
        

end
