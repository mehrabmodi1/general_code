clear all
close all

dataset_list_paths = [...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONAlpha1.xls'};...
                    ];
            
suppress_plots = 0;
[del, odor_names1] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
[del, ]
global color_vec;                
a = colormap('bone');
global greymap
greymap = flipud(a);
fly_n = 0;
script_name = mfilename;

for list_n = 1:size(dataset_list_paths, 1)
    curr_dir_list_path = dataset_list_paths{list_n, 1};
    [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
    n_dirs = size(dir_list, 1);
    dataset_list_name = findstr(curr_dir_list_path, 'El_');
    dataset_list_name = curr_dir_list_path((dataset_list_name + 1):(end - 4));
    
    dataset_list_name(1) = [];
    flies_resp_size_mat = [];
    saved_resp_sizes_all = [];
    
    %loop to go through all experiment datasets listed in list file
    for dir_n = 1:n_dirs
       fly_n = fly_n + 1;
              
       saved_an_results.scriptname = mfilename('fullpath');
       curr_dir = [dir_list{dir_n, 1}, '\'];
       curr_dir = manage_base_paths(curr_dir, 2);
       
       tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
       tif_times = tif_times.time_stamps;
       [stim_mat, stim_mat_simple, column_heads, color_vec] = load_params_trains(curr_dir, tif_times);    %reading in trial stimulus parameters after matching time stamps to F traces
       odor_colors = [color_vec(3, :); color_vec(3, :).*0.75; color_vec(2, :)];
       
       %Reading in experimental parameters
        %odor_list = unique(stim_mat_simple(:, 2) );
        odor_list = [3, 10, 11];
        n_odors = length(odor_list);
        odor_dur_list = unique(stim_mat_simple(:, 3) );
        n_od_durs = length(odor_dur_list);
        n_trains = max(stim_mat_simple(:, 11));
        saved_an_results.odor_list = odor_list;
        saved_an_results.odor_dur_list = odor_dur_list;
        fly_resp_size_vec = zeros(1, length(odor_list).*2) + nan;
        
        cd(curr_dir);
        tif_name = dir('*.tif');
        stack_obj = ScanImageTiffReader([curr_dir, tif_name(1).name]);
        [frame_time, zoom, n_chans, PMT_offsets] = SI_tif_info(stack_obj);
        
        
        %loading extracted raw fluorescence data matrices written by raw_dff_extractor
        raw_data_mat = load([curr_dir 'extracted_raw_data_mat.mat']);
        raw_data_mat = raw_data_mat.raw_data_mat;           %raw F traces extracted from ROIs
        raw_data_mat_orig = raw_data_mat;
        raw_data_mat = raw_data_mat(:, :, stim_mat_simple(:, 1));       %making sure only time-stamp matched trials are used for further analysis
        n_cells = size(raw_data_mat, 2);
        
        %calculating dF/F traces from raw data
        filt_time = 0.5;            %in s, the time window for boxcar filter for generating filtered traces
        [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, curr_dir);
        del = find(dff_data_mat_f < -1);
        dff_data_mat_f(del) = -1;       %forcing crazy values to sane ones
        