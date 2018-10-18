clear all
close all

%root_path = 'C:\Data\Data\Analysed_data\Analysis_results\Yoshi_PaBaEl\';
root_path = 'C:\Data\Data\Analysed_data\Analysis_results\Yoshi_PaBaEl_lowflux\';

cd(root_path);
result_file_list = dir('*.mat');
[del, odor_names] = xlsread('C:\Data\Code\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);

pair_names = [{'PA-BA'}, {'PA-EL'}, {'BA-EL'}];
odor_names = [{'PA'}, {'BA'}, {'EL'}];

a = colormap('bone');
global greymap
greymap = flipud(a);
color_vec = setup_std_color_vec(3);

suppress_plots = 1;

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
    int_resps_mat_1s = saved_an_results.int_resps_1s;
    non_int_resps_mat_1s = saved_an_results.non_int_resps_1s;
    int_resps_mat_60s = saved_an_results.int_resps_60s;
    non_int_resps_mat_60s = saved_an_results.non_int_resps_60s;
    stim_frs = saved_an_results.stim_frs_saved;
    
    for dur_n = 1:2
        mean_color = [0.84, 0.38, 0.3];
        mat = squeeze(sparseness_mat(:, dur_n, :));
        figure(1)
        scattered_dot_plot(mat', 1, 0.5, 0.5, 8, [0.3, 0.3, 0.3], 0, [], [1, 1, 1], odor_names, 1, mean_color);
        ylabel('responder fraction')
        title([dataset_name, ', stimulus duration ' num2str(od_dur_list(dur_n)), 's'])
        axis([0, 3, 0, 0.5]);
        fig_wrapup(1)
                
        %plotting intersections
        mat = squeeze(intersection_mat(:, dur_n, :));
        mat(isnan(mat)) = 0;
        figure(2)
        scattered_dot_plot(mat', 2, 0.5, 0.5, 8, [0.3, 0.3, 0.3], 1, [], [0.5, 0.5, 0.5], pair_names, 1, mean_color);
        ylabel('repr. overlap (fold expected)')
        title([dataset_name, ', stimulus duration ' num2str(od_dur_list(dur_n)), 's'])
        axis([0, 3, 0, 15]);
        fig_wrapup(2)
        
        %plotting intersections norm to n_cells
        mat = squeeze(intersection_mat_n(:, dur_n, :));
        mat(isnan(mat)) = 0;
        figure(3)
        scattered_dot_plot(mat', 3, 0.5, 0.5, 8, [0.5, 0.5, 0.5], 1, [], [0.5, 0.5, 0.5], pair_names, 1, mean_color);
        ylabel('repr. overlap (frac. cells)')
        title([dataset_name, ', stimulus duration ' num2str(od_dur_list(dur_n)), 's'])
        axis([0, 3, 0, 0.5]);
        fig_wrapup(3)
        
        %plotting non-overlap norm to n_cells
        mat = squeeze(non_int_mat(:, dur_n, :));
        mat(isnan(mat)) = 0;
        figure(4)
        scattered_dot_plot(mat', 4, 0.5, 0.5, 8, [0.5, 0.5, 0.5], 1, [], [0.5, 0.5, 0.5], pair_names, 1, mean_color);
        ylabel('repr. nonoverlap (frac. cells)')
        title([dataset_name, ', stimulus duration ' num2str(od_dur_list(dur_n)), 's'])
        axis([0, 3, 0, 0.5]);
        fig_wrapup(4)
        
        if suppress_plots == 0
            keyboard
        else
        end
    
    end
    
    %plotting distributions of significant response sizes
    resp_vec = saved_an_results.pk_responses;
    bins = 0:0.1:4; 
    [counts] = hist(resp_vec, bins);
    counts = counts./sum(counts);
    figure(5)
    if isempty(findstr(dataset_name, 'AlphaBeta')) == 0
        line_col = [0.9, 0.4, 0.3];
        resp_vec1 = resp_vec;
    elseif isempty(findstr(dataset_name,'Gamma')) == 0
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
    
    %plotting resp mats of cells sorted by time of center of mass
        
    %sorting cells by times of centers of mass of response vectors
    int_1s_CoM_vec = sum(int_resps_mat_1s(stim_frs(1, 1):(stim_frs(1, 2) + 200), :, 1), 1, 'omitnan');
    int_resps_mat_1s_1 = sort_by_vec(int_resps_mat_1s(:, :, 1), int_1s_CoM_vec, 1);
    int_resps_mat_1s_2 = sort_by_vec(int_resps_mat_1s(:, :, 2), int_1s_CoM_vec, 1);
    int_60s_CoM_vec = sum(int_resps_mat_60s(stim_frs(2, 1):(stim_frs(2, 2) + 50), :, 1), 1, 'omitnan');
    int_resps_mat_60s_1 = sort_by_vec(int_resps_mat_60s(:, :, 1), int_60s_CoM_vec, 1);
    int_resps_mat_60s_2 = sort_by_vec(int_resps_mat_60s(:, :, 2), int_60s_CoM_vec, 1);
    
    non_int_1s_CoM_vec = sum(non_int_resps_mat_1s(stim_frs(1, 1):(stim_frs(1, 2) + 200), :, 1), 1, 'omitnan');
    non_int_resps_mat_1s_1 = sort_by_vec(non_int_resps_mat_1s(:, :, 1), non_int_1s_CoM_vec, 1);
    non_int_resps_mat_1s_2 = sort_by_vec(non_int_resps_mat_1s(:, :, 2), non_int_1s_CoM_vec, 1);
    
    non_int_60s_CoM_vec = sum(non_int_resps_mat_60s(stim_frs(2, 1):(stim_frs(2, 2) + 50), :, 1), 1, 'omitnan');
    non_int_resps_mat_60s_1 = sort_by_vec(non_int_resps_mat_60s(:, :, 1), non_int_60s_CoM_vec, 1);
    non_int_resps_mat_60s_2 = sort_by_vec(non_int_resps_mat_60s(:, :, 2), non_int_60s_CoM_vec, 1);
    
    figure(6)
    imagesc(int_resps_mat_1s_1', [0, 2])
    colormap(greymap)
    ylabel('cell number')
    set_xlabels_time(6, 0.099, 2)
    fig_wrapup(6)
    add_stim_bar(6, stim_frs(1, :), [0.5, 0.5, 0.5]);
    
    
    figure(7)
    imagesc(int_resps_mat_1s_2', [0, 2])
    colormap(greymap)
    ylabel('cell number')
    set_xlabels_time(7, 0.099, 2)
    fig_wrapup(7)
    add_stim_bar(7, stim_frs(1, :), [0.5, 0.5, 0.5]);
    
    
    figure(8)
    imagesc(int_resps_mat_60s_1', [0, 2])
    colormap(greymap)
    ylabel('cell number')
    set_xlabels_time(8, 0.099, 2)
    fig_wrapup(8)
    add_stim_bar(8, stim_frs(2, :), [0.5, 0.5, 0.5]);
    
    
    figure(9)
    imagesc(int_resps_mat_60s_2', [0, 2])
    colormap(greymap)
    ylabel('cell number')
    set_xlabels_time(9, 0.099, 2)
    %cb = colorbar('eastoutside');
    fig_wrapup(9)
    add_stim_bar(9, stim_frs(2, :), [0.5, 0.5, 0.5]);
    
    figure(10)
    imagesc(non_int_resps_mat_1s_1', [0, 2])
    colormap(greymap)
    ylabel('cell number')
    set_xlabels_time(10, 0.099, 2)
    fig_wrapup(10)
    add_stim_bar(10, stim_frs(1, :), [0.5, 0.5, 0.5]);
    
    
    figure(11)
    imagesc(non_int_resps_mat_1s_2', [0, 2])
    colormap(greymap)
    ylabel('cell number')
    set_xlabels_time(11, 0.099, 2)
    fig_wrapup(11)
    add_stim_bar(11, stim_frs(1, :), [0.5, 0.5, 0.5]);
    
    
    figure(12)
    imagesc(non_int_resps_mat_60s_1', [0, 2])
    colormap(greymap)
    ylabel('cell number')
    set_xlabels_time(12, 0.099, 2)
    fig_wrapup(12)
    add_stim_bar(12, stim_frs(2, :), [0.5, 0.5, 0.5]);
    
    
    figure(13)
    imagesc(non_int_resps_mat_60s_2', [0, 2])
    colormap(greymap)
    ylabel('cell number')
    set_xlabels_time(13, 0.099, 2)
    %cb = colorbar('eastoutside');
    fig_wrapup(13)
    add_stim_bar(13, stim_frs(2, :), [0.5, 0.5, 0.5]);
    
     
    close figure 6
    close figure 7
    close figure 8
    close figure 9
    close figure 10
    close figure 11
    close figure 12
    close figure 13
end
[h, p] = ttest2(resp_vec1, resp_vec2)   %testing for sig difference in ave resp size