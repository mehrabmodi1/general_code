clear all
close all

%PICK UP THREAD HERE
%Figure out how to plot saved int and non_int traces

root_path = 'C:\Data\Data\Analysed_data\Analysis_results\Yoshi_PaBaEl\';
cd(root_path);
result_file_list = dir('*.mat');
[del, odor_names] = xlsread('C:\Data\Code\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);

pair_names = [{'PA-BA'}, {'PA-EL'}, {'BA-EL'}];
odor_names = [{'PA'}, {'BA'}, {'EL'}];

for result_file_n = 1:size(result_file_list, 1)
    saved_an_results = load([root_path, result_file_list(result_file_n).name]);
    saved_an_results = saved_an_results.saved_an_results;
    odor_list = saved_an_results.odor_list;    
    od_dur_list = saved_an_results.odor_dur_list;
    dataset_namei = findstr(result_file_list(result_file_n).name, '_an');
    dataset_name = result_file_list(result_file_n).name(1:(dataset_namei - 1));
    
    %plotting sparsenesses
    sparseness_mat = saved_an_results.sparsenesses;
    intersection_mat = saved_an_results.sig_intersections;
    intersection_mat_n = saved_an_results.sig_intersections_n;
    non_int_mat = saved_an_results.sig_nonintersections;
    
    for dur_n = 1:2
        mean_color = [0.84, 0.38, 0.3];
        mat = squeeze(sparseness_mat(:, dur_n, :));
        figure(1)
        scattered_dot_plot(mat', 1, 0.5, 0.5, 8, [0.3, 0.3, 0.3], 0, [1, 1, 1], odor_names, 1, mean_color);
        ylabel('responder fraction')
        title([dataset_name, ', stimulus duration ' num2str(od_dur_list(dur_n)), 's'])
        axis([0, 3, 0, 0.5]);
        fig_wrapup(1)
                
        %plotting intersections
        mat = squeeze(intersection_mat(:, dur_n, :));
        mat(isnan(mat)) = 0;
        figure(2)
        scattered_dot_plot(mat', 2, 0.5, 0.5, 8, [0.3, 0.3, 0.3], 1, [0.5, 0.5, 0.5], pair_names, 1, mean_color);
        ylabel('repr. overlap (fold expected)')
        title([dataset_name, ', stimulus duration ' num2str(od_dur_list(dur_n)), 's'])
        axis([0, 3, 0, 15]);
        fig_wrapup(2)
        
        %plotting intersections norm to n_cells
        mat = squeeze(intersection_mat_n(:, dur_n, :));
        mat(isnan(mat)) = 0;
        figure(3)
        scattered_dot_plot(mat', 3, 0.5, 0.5, 8, [0.5, 0.5, 0.5], 1, [0.5, 0.5, 0.5], pair_names, 1, mean_color);
        ylabel('repr. overlap (frac. cells)')
        title([dataset_name, ', stimulus duration ' num2str(od_dur_list(dur_n)), 's'])
        axis([0, 3, 0, 0.5]);
        fig_wrapup(3)
        
        %plotting non-overlap norm to n_cells
        mat = squeeze(non_int_mat(:, dur_n, :));
        mat(isnan(mat)) = 0;
        figure(4)
        scattered_dot_plot(mat', 4, 0.5, 0.5, 8, [0.5, 0.5, 0.5], 1, [0.5, 0.5, 0.5], pair_names, 1, mean_color);
        ylabel('repr. nonoverlap (frac. cells)')
        title([dataset_name, ', stimulus duration ' num2str(od_dur_list(dur_n)), 's'])
        axis([0, 3, 0, 0.5]);
        fig_wrapup(4)
        
    end
    
    %plotting distributions of significant response sizes
    resp_vec = saved_an_results.pk_responses;
    bins = 0:0.1:4; 
    [counts] = hist(resp_vec, bins);
    counts = counts./sum(counts);
    figure(5)
    if strcmp(dataset_name, 'AlphaBeta') == 1
        line_col = [0.9, 0.4, 0.3];
        resp_vec1 = resp_vec;
    elseif strcmp(dataset_name,'Gamma') == 1
        line_col = [0.3, 0.9, 0.4];
        resp_vec2 = resp_vec;
    else
    end
    
    plot(bins, counts, 'Color', line_col, 'lineWidth', 3)
    hold on
    legend('Alpha/Beta', 'Gamma', 'Location', 'northeast')
    xlabel('mean response size (dF/F)')
    ylabel('frac. significant responses')
    fig_wrapup(5)
    
end
[h, p] = ttest2(resp_vec1, resp_vec2)