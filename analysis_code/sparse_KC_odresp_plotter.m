clear all
close all

dataset_list_paths = [...
                      
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\KC_sparse_odor_resps_set1.xls'};...

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
    
    saved_responses = [];           %8 cols per row, one row per axon; 
    saved_responses_mean = [];
    
    %loop to go through all experiment datasets listed in list file
    for dir_n = 1:n_dirs
        fly_n = fly_n + 1;
              
        saved_an_results.scriptname = mfilename('fullpath');
        curr_dir = [dir_list{dir_n, 1}, '\'];
        curr_dir = manage_base_paths(curr_dir, 2);
       
        tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
        tif_times = tif_times.time_stamps;
        [stim_mat, stim_mat_simple, column_heads, color_vec] = load_params_trains_modular(curr_dir, tif_times);    %reading in trial stimulus parameters after matching time stamps to F traces
                
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
        filt_time = 5;            %in s, the time window for boxcar filter for generating filtered traces
        [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, curr_dir);
%         del = find(dff_data_mat_f < -1);
%         dff_data_mat_f(del) = -1;       %forcing crazy values to sane ones
                
       
        
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
        ROI_reinforced_vec = [];
        for ROI_n = 2:size(ROI_mat, 3)
            sum_mat = ROI_mat(:, :, ROI_n) + ROI_mat(:, :, 1);
            if max(max(sum_mat)) == 2
                ROI_reinforced_vec = [ROI_reinforced_vec; ROI_n];           %ROIs that received DAN reinforcement
            else
            end
                
        end
           
        
        %manually defining sets of ROIs (to label ROI sets that lie along an axon)
        %[ROI_sets, marked_coords] = define_ROI_sets(ROI_mat, 1);
        if exist([curr_dir, 'ROI_sets.mat']) == 2
            ROI_sets = load([curr_dir, 'ROI_sets.mat']);
            ROI_sets = ROI_sets.ROI_sets;
        else
            [ROI_sets] = define_roi_groups(ROI_mat, ROI_mat(:, :, 1), 1);
            save([curr_dir, 'ROI_sets.mat'], 'ROI_sets');
        end
        
        %plotting traces and building matrix of response sizes
        n_axons = size(ROI_sets, 1);
        curr_od_list = odor_list_olf1;
        
        
        for axon_n = 0:(n_axons - 1)
            resp_vec = [];
            resp_vec_mean = [];
            
            curr_ROIs = ROI_sets{(axon_n + 1), 1};
            %making sure paired ROI (ROI in DAN axons) is first
            paired_ROI = intersect(curr_ROIs, ROI_reinforced_vec);
            curr_ROIs(curr_ROIs == paired_ROI) = [];
            curr_ROIs = [paired_ROI, curr_ROIs];

            for ROI_ni = 1:min([length(curr_ROIs), 2])
                ROI_n = curr_ROIs(ROI_ni);
                for odor_ni = 1:size(curr_od_list, 1)
                    odor_n = curr_od_list(odor_ni);
                    od_name = odor_names1{odor_n, 1};
                    curr_trs = find(stim_mat_simple(:, od_olf1_col_n) == odor_n);
                    stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
                    stim_frs = stim_frs{1};
                    
                    curr_traces = squeeze(dff_data_mat_f(:, ROI_n, curr_trs));
                    mean_trace = mean(curr_traces, 2, 'omitnan');
                    se_trace = std(curr_traces, [], 2, 'omitnan')./sqrt(size(curr_traces, 2));
                    max_resp = max(mean_trace(stim_frs(1):(stim_frs(2) + round(7./frame_time))) );     %taking the max in an extended window allows for on/sus/off responses to be captured. trace is filtered.
                    mean_resp = mean(mean_trace(stim_frs(1):(stim_frs(2) + round(7./frame_time))) );     %taking the max in an extended window allows for on/sus/off responses to be captured. trace is filtered.
                             
                    resp_vec = [resp_vec, max_resp];
                    resp_vec_mean = [resp_vec_mean, mean_resp];
                    
                    
                    if suppress_plots == 0
                        %plotting traces
                        fig_n = 1;
                        figure(fig_n)
                        plot(curr_traces, 'lineWidth', 2, 'Color', color_vec(odor_ni, :))
                        
                        set_xlabels_time(fig_n, frame_time, 25);
                        ylabel(['ROI ', int2str(ROI_n), ' ', od_name, ' resp. (dF/F)']);
                        fig_wrapup(fig_n, script_name);
                        add_stim_bar(fig_n, stim_frs, [0.75, 0.75, 0.75]);
                        
                        fig_n = fig_n + 1;
                        figure(fig_n)
                        shadedErrorBar([], mean_trace, se_trace, {'Color', color_vec(odor_ni, :)}, 1);
                        ylabel(['ROI ', int2str(ROI_n), ' ', od_name, ' resp. (dF/F)']);
                        set_xlabels_time(fig_n, frame_time, 25);
                        fig_wrapup(fig_n, script_name);
                        add_stim_bar(fig_n, stim_frs, [0.75, 0.75, 0.75]);
                        
                        keyboard
                        
                    else
                    end
                    
                    close figure 1
                    close figure 2
                end

                

                
                

            end
            resp_vec_pad = nan(1, 8);
            resp_vec_pad(1, 1:length(resp_vec)) = resp_vec;
            resp_vec = resp_vec_pad;
            resp_vec_pad = nan(1, 8);
            resp_vec_pad(1, 1:length(resp_vec_mean)) = resp_vec_mean;
            resp_vec_mean = resp_vec_pad;
            saved_responses = [saved_responses; resp_vec];
            saved_responses_mean = [saved_responses_mean; resp_vec_mean];
        end
        
    end
    keyboard
end