clear all
close all

dataset_list_paths = [%{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_d5HT1b_Gamma.xls'};...
                      {'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_d5HT1b_AlphaBeta.xls'} ...
                      ];

for list_n = 1:size(dataset_list_paths, 1)
    curr_dir_list_path = dataset_list_paths{list_n, 1};
    [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
    n_dirs = size(dir_list, 1);
   
    %loop to go through all experiment datasets listed in list file
    for dir_n = 1:n_dirs
        curr_dir = [dir_list{dir_n, 1}, '\'];