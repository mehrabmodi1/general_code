clear all
close all

dataset_list_paths = [...
                      
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_KC_sparse_plasticity_set1.xls'};...

                      ];
            
suppress_plots = 0;
[del, odor_names1] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
[del, odor_names2] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList_olf2.xls', 1);

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
    dataset_list_name = findstr(curr_dir_list_path, 'list_');
    dataset_list_name = curr_dir_list_path((dataset_list_name + 5):(end - 4));
    
    saved_responses = zeros(8, n_dirs) + nan;     %2 odours, 2 ROIs, 2 conditions - pre and post
    
    %loop to go through all experiment datasets listed in list file
    for dir_n = 1:n_dirs
        fly_n = fly_n + 1;
              
        saved_an_results.scriptname = mfilename('fullpath');
        curr_dir = [dir_list{dir_n, 1}, '\'];
        curr_dir = manage_base_paths(curr_dir, 2);
       
        tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
        tif_times = tif_times.time_stamps;
        [stim_mat, stim_mat_simple, column_heads, color_vec] = load_params_trains_modular(curr_dir, tif_times);    %reading in trial stimulus parameters after matching time stamps to F traces
        fore_colour = color_vec(1, :);
        distr_colour = color_vec(2, :);
        ctrl_colour = color_vec(3, :);
        
        %identifying stim_mat_simple col numbers
        led_on_col_n = find_stim_mat_simple_col('led_on', column_heads);            %identifying relevant column number in stim_mat_simple
        od_olf1_col_n = find_stim_mat_simple_col('odor_n', column_heads);           %identifying relevant column number in stim_mat_simple
        od_olf2_col_n = find_stim_mat_simple_col('odour_olf2', column_heads);       %identifying relevant column number in stim_mat_simple
        dur_olf1_col_n = find_stim_mat_simple_col('duration', column_heads);        %identifying relevant column number in stim_mat_simple
        dur_olf2_col_n = find_stim_mat_simple_col('duration_olf2', column_heads);   %identifying relevant column number in stim_mat_simple
        od_col_ns = [od_olf1_col_n, od_olf2_col_n];
        dur_col_ns = [dur_olf1_col_n, dur_olf2_col_n];        
        
        %Reading in experimental parameters
        odor_list_olf1 = unique(stim_mat_simple(:, od_olf1_col_n) );
        n_odors_olf1 = length(odor_list_olf1);
        odor_dur_list_olf1 = unique(stim_mat_simple(:, dur_olf1_col_n) );
        n_od_durs_olf1 = length(odor_dur_list_olf1);
        
        cd(curr_dir);
        tif_name = dir('*.tif');
        stack_obj = ScanImageTiffReader([curr_dir, tif_name(1).name]);
        [frame_time, zoom, n_chans, PMT_offsets] = SI_tif_info(stack_obj);
        
        
        %loading extracted raw fluorescence data matrices written by raw_dff_extractor
        raw_data_mat = load([curr_dir 'extracted_raw_data_mat.mat']);
        raw_data_mat = raw_data_mat.raw_data_mat;           %raw F traces extracted from ROIs
        raw_data_mat_orig = raw_data_mat;
        ROI_mat = load([curr_dir, '\ROI_mat.mat']);
        ROI_mat = ROI_mat.ROI_mat;
        tif_n_col_n = find_stim_mat_simple_col('matched_tif_n', column_heads);
        raw_data_mat = raw_data_mat(:, :, stim_mat_simple(:, tif_n_col_n));       %making sure only time-stamp matched trials are used for further analysis
        n_cells = size(raw_data_mat, 2);
        
        %calculating dF/F traces from raw data
        filt_time = 0.5;            %in s, the time window for boxcar filter for generating filtered traces
        [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, curr_dir);
        del = find(dff_data_mat_f < -1);
        dff_data_mat_f(del) = -1;       %forcing crazy values to sane ones
                
        pairing_tr_n = find(stim_mat_simple(:, led_on_col_n) == 1);
        paired_odor_n = stim_mat_simple(pairing_tr_n, od_olf1_col_n);
        unpaired_odor_n = odor_list_olf1(odor_list_olf1~=paired_odor_n);
        
        %identifying pre and post pairing trial sets
        pairing_tr = find(stim_mat_simple(:, 18) == 1);
        paired_odor_vec(dir_n, 1) = stim_mat(pairing_tr).odours;     %keeping track of which odor was the paired one        
        
        %identifying the ROI that marks out the region of DAN innervation
        n_ROIs = size(ROI_mat, 3);
        max_area = [0, 0];
        for ROI_n = 1:n_ROIs
            curr_area = sum(sum(ROI_mat(:, :, ROI_n)));
            if curr_area > max_area(1, 1)
                max_area = [curr_area, ROI_n];
            else
            end
        end
        big_ROI = ROI_mat(:, :, max_area(1, 2));
        ROI_mat(:, :, max_area(1, 2)) = [];
        ROI_mat = cat(3, big_ROI, ROI_mat);
        
                
        %identifying ROIs within DAN innervation zone
        ROI_reinforced_vec = zeros(n_ROIs, 1);
        for ROI_n = 2:n_ROIs
            sum_mat = ROI_mat(:, :, ROI_n) + ROI_mat(:, :, 1);
            if max(max(sum_mat)) < 2
                ROI_reinforced_vec(ROI_n, 1) = 0;           %ROIs that did not receive DAN reinforcement
            elseif max(max(sum_mat)) == 2           
                ROI_reinforced_vec(ROI_n, 1) = 1;           %ROIs that received DAN reinforcement
            else
            end
                
        end
                
        
        %manually defining sets of ROIs (to label ROI sets that lie along an axon)
        %[ROI_sets, marked_coords] = define_ROI_sets(ROI_mat, 1);
        [ROI_sets] = define_roi_groups(ROI_mat, ROI_mat(:, :, 1));
        keyboard
        
        n_axons = size(ROI_sets, 1);
        curr_odor = odor_list_olf1();
            for axon_n = 0:(n_axons - 1)
                curr_ROIs = ROI_sets{axon_n, 1};

                for ROI_ni = 1:size(curr_ROIs)
                    ROI_n = curr_ROIs(ROI_ni);

                    %curr_traces = dff_data_mat_f(:, );

                end



                fig_n = (axon_n + 1)*3 + 1;
                figure(fig_n)


            end
        
   
    end
    
   
end