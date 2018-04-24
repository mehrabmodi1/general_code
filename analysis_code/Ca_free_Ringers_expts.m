clear all
close all

direc_list_path = 'C:\Data\Data\Raw_data\dataset_lists\dataset_list_Cafree_Ringers_expt.xls';
[del, direc_list] = xlsread(direc_list_path, 1);
n_direcs = size(cell_direc_list, 1);

for direc_n = 1:n_direcs
    curr_dir = direc_list{direc_n, 1};
    
    dir_contents = dir_date_sorted(curr_dir, '*.tif');
    
    
    keyboard
    
end