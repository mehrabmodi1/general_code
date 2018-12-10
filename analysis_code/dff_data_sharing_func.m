clear all
close all

dataset_list_paths = [%{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_d5HT1b_Gamma.xls'};...
                      %{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_c739_AlphaBeta.xls'} ...
                      %{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_d5HT1b_Gamma_low_conc.xls'};...
                      %{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_c739_AlphaBeta_low_conc.xls'};...
                      %{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_c305a_AlphapBetap_low_conc.xls'};...
                      {'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_d5HT1b_Gamma_set2.xls'};...
                      {'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_c739_AlphaBeta_set2.xls'}...
                      ];

suppress_plots = 1;
[del, odor_names] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
a = colormap('bone');
global greymap
greymap = flipud(a);


for list_n = 1:size(dataset_list_paths, 1)
    curr_dir_list_path = dataset_list_paths{list_n, 1};
    [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
    n_dirs = size(dir_list, 1);
    dataset_list_name = findstr(curr_dir_list_path, 'El_');
    dataset_list_name = curr_dir_list_path((dataset_list_name + 1):(end - 4));
    dataset_list_name(1) = [];
    
    
    for dir_n = 1:n_dirs
        saved_an_results.scriptname = mfilename('fullpath');
        curr_dir = [dir_list{dir_n, 1}, '\'];
        tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
        tif_times = tif_times.time_stamps;
        [stim_mat, stim_mat_simple, column_heads, color_vec, g_tr_list] = load_params_trains(curr_dir, tif_times);    %reading in trial stimulus parameters after matching time stamps to F traces
                
        %Reading in experimental parameters
        odor_list = unique(stim_mat_simple(:, 2) );
        n_odors = length(odor_list);
        odor_dur_list = unique(stim_mat_simple(:, 3) );
        n_od_durs = length(odor_dur_list);
        n_trains = max(stim_mat_simple(:, 11));
        saved_an_results.odor_list = odor_list;
        saved_an_results.odor_dur_list = odor_dur_list;
        
        cd(curr_dir);
        tif_name = dir('*.tif');
        stack_obj = ScanImageTiffReader([curr_dir, tif_name(1).name]);
        [frame_time, zoom, n_chans, PMT_offsets] = SI_tif_info(stack_obj);

       
        %loading extracted raw fluorescence data matrices written by raw_dff_extractor
        raw_data_mat = load([curr_dir 'extracted_raw_data_mat.mat']);
        raw_data_mat = raw_data_mat.raw_data_mat;           %raw F traces extracted from ROIs
        raw_data_mat_orig = raw_data_mat;
        raw_data_mat = raw_data_mat(:, :, stim_mat_simple(:, 1));       %making sure only time-stamp matched trials are used for further analysis
        bad_tr_list = 1:1:size(raw_data_mat, 3);
        bad_tr_list(g_tr_list) = [];
        raw_data_mat(:, :, bad_tr_list) = nan;
        n_cells = size(raw_data_mat, 2);
        
        
        %calculating dF/F traces from raw data
        filt_time = 1;            %in s, the time window for boxcar filter for generating filtered traces
        [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, curr_dir);
        del = find(dff_data_mat_f < -1);
        dff_data_mat_f(del) = -1;       %forcing crazy values to sane ones
        %identifying significantly responsive cells
        [resp_areas, sig_trace_mat, sig_trace_mat_old, sig_cell_mat] = cal_sig_responses_res(dff_data_mat, stim_mat, stim_mat_simple, curr_dir, frame_time);

        %Running data quality control checks
%         sig_cell_mat_old = sig_cell_mat;
%         [sig_cell_mat, all_bad_trs] = cell_data_quality_control(dff_data_mat_f, stim_mat, stim_mat_simple, sig_cell_mat, 1, frame_time);
%         dff_data_mat(:, :, all_bad_trs) = nan;
        %disp([num2str((length(all_bad_trs)./size(dff_data_mat, 3)).*100) ' percent of trials were auto-identified as bad and removed.']);
        sig_cells = find(sum(sum(sig_cell_mat, 3), 2) > 0);         %list of all cells significant for any odor for any duration
       
        %visualising odor sorted responses
        [resp_areas_sorted, fig_h, c_vals_mat] = data_quality_resp_matrix(resp_areas, stim_mat_simple, sig_cell_mat, 1, 1);
        pop_c_vals(dir_n, :) = c_vals_mat;
        
        
        %plotting ave KC pop resp traces to each odor
        KCpop_odor_resp_plotter(dff_data_mat_f, stim_mat_simple, stim_mat, sig_cell_mat, frame_time, 1);
        keyboard
        close all
        %saving stuff for sharing data
        if list_n == 1
            share_path = ['C:\Data\Data\Analysed_data\data_sharing\Gamma_set2\fly' int2str(dir_n)];
        elseif list_n == 2
            share_path = ['C:\Data\Data\Analysed_data\data_sharing\AlphaBeta_set2\fly' int2str(dir_n)];
        else
        end
        mkdir(share_path);
        save([share_path, '\dFF_data.mat'], 'dff_data_mat_f');
        save([share_path, '\sig_cells.mat'], 'sig_cell_mat');
        save([share_path, '\stim_mat.mat'], 'stim_mat')
        
    end
    keyboard
end