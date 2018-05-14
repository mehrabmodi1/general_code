clear all
close all

[del, odor_names] = xlsread('C:\Data\Code\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
dir_list_path = 'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_odpulse_freq_series.xls';

[del, dir_list] = xlsread(dir_list_path, 1);        %list of Suite2P results directories
n_dirs = size(dir_list, 1);

%loop to go through all experiment datasets listed in list file
for dir_n = 1:n_dirs
    curr_dir = [dir_list{dir_n, 1}, '\'];
    tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
    tif_times = tif_times.time_stamps;
    [stim_mat, stim_mat_simple, column_heads, color_vec] = load_params_trains(curr_dir, tif_times);    %reading in trial stimulus parameters after matching time stamps to F traces
    
    %Reading in experimental parameters
    odor_list = unique(stim_mat_simple(:, 2) );
    n_odors = length(odor_list);
    odor_dur_list = unique(stim_mat_simple(:, 3) );
    n_od_durs = length(odor_dur_list);
    n_trains = max(stim_mat_simple(:, 11));
    
    cd(curr_dir);
    tif_name = dir('*.tif');
    stack_obj = ScanImageTiffReader([curr_dir, tif_name(1).name]);
    
    keyboard
    
    
    
end