clear all
close all

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
    for dur_n = 1:2
        mat = squeeze(sparseness_mat(:, dur_n, :));
        figure(1)
        scattered_dot_plot(mat', 1, 0.5, 0.5, 8, [0.3, 0.3, 0.3], 0, [1, 1, 1], odor_names);
        ylabel('responder fraction')
        title([dataset_name, ', stimulus duration ' num2str(od_dur_list(dur_n)), 's'])
        axis([0, 3, 0, 1]);
        fig_wrapup(1)
        %keyboard
    end
    
    %plotting intersections
    intersection_mat = saved_an_results.sig_intersections;
    for dur_n = 1:2
        mat = squeeze(intersection_mat(:, dur_n, :));
        mat(isnan(mat)) = 0;
        figure(1)
        scattered_dot_plot(mat', 1, 0.5, 0.5, 8, [0.3, 0.3, 0.3], 1, [0.5, 0.5, 0.5], pair_names);
        ylabel('repr. overlap (fold expected)')
        title([dataset_name, ', stimulus duration ' num2str(od_dur_list(dur_n)), 's'])
        axis([0, 3, 0, 15]);
        fig_wrapup(1)
        keyboard
    end
    
    
   keyboard 
end