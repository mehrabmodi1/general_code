clear all
close all

dataset_list_paths = [%{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_d5HT1b_Gamma.xls'};...
                      %{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_c739_AlphaBeta.xls'} ...
                      %{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_d5HT1b_Gamma_low_conc.xls'};...
                      %{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_c739_AlphaBeta_low_conc.xls'};...
                      %{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_c305a_AlphapBetap_low_conc.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_Yoshi_PaBaEl_d5HT1b_Gamma_set2.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_Yoshi_PaBaEl_c739_AlphaBeta_set2.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_rand_trains.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_KC_Ca_alpha1T_set2.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_KC_c739_PABAEL_201908set.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_KC_d5HT1b_PABAEL_201908set.xls'}...                      
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_MBONG2_PaBaEl_handover_starved_halfAra_prehabituated_strongUS.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_c739KC_PaBaEl_handover_prehabituated_trimmed.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\c305aKC_PaBaEl_handover_rtrain_prehabituated.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\KCd5HT1b_Gamma_PaBaEl_handover_rtrain_prehabituated.xls'};...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\c739KC_PaBaEl_handover_rtrain_prehabituated.xls'};...
                      
                      ];

%saving stuff for sharing data
share_path_base = 'C:\Data\Data\Analysed_data\data_sharing\param_fitting\';
dataset_list_name = 'ABKC_hover_rtrain_prehabituated\';    


dataset_type = 3; %2 - manualROI, 3 - Suite2P
[del, odor_names] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
a = colormap('bone');
global greymap
greymap = flipud(a);
script_name = mfilename;

for list_n = 1:size(dataset_list_paths, 1)
    curr_dir_list_path = dataset_list_paths{list_n, 1};
    [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
    n_dirs = size(dir_list, 1);
%     dataset_list_name = findstr(curr_dir_list_path, 'dataset_list');
%     dataset_list_name = curr_dir_list_path((dataset_list_name + 13):(end - 4));
   
   
    
    for dir_n = 1:n_dirs
        saved_an_results.scriptname = mfilename('fullpath');
        curr_dir = [dir_list{dir_n, 1}, '\'];
        curr_dir = manage_base_paths(curr_dir, dataset_type);
        curr_dir = raw_direc_with_1(curr_dir);
        curr_dir = [curr_dir, '\'];
        tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
        tif_times = tif_times.time_stamps;
        [stim_mat, stim_mat_simple, column_heads, color_vec, good_tr_list, params_orig] = load_params_trains_modular(curr_dir, tif_times);    %reading in trial stimulus parameters after matching time stamps to F traces%[stim_mat, stim_mat_simple, column_heads, color_vec, g_tr_list] = load_params_trains(curr_dir, tif_times);    %reading in trial stimulus parameters after matching time stamps to F traces
        
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
        tif_n_col_n = find_stim_mat_simple_col('matched_tif_n', column_heads);
        
        %inserting dummy trials (all nans) into raw_data_mat for pairing trials for
        %which no corress .tifs were acquired
        raw_data_mat = match_up_rawmat_matchedtrs(raw_data_mat, stim_mat_simple, tif_n_col_n);
        
        bad_tr_list = 1:1:size(raw_data_mat, 3);
        bad_tr_list(good_tr_list) = [];
        raw_data_mat(:, :, bad_tr_list) = nan;
        n_cells = size(raw_data_mat, 2);
        
        %creating a list of trial numbers to keep track of which trials
        %were included for KC param fitting.
        tr_list = 1:size(raw_data_mat, 3);
        
        %since this is for KC response parameter extraction, only saving
        %the longest duration, single-odor, single-pulse trials for all odors. Also
        %combining olf1 and olf2 paramters as follows: olf2 od_ns will be
        %listed as olf1_od_ns starting from od_n 13 and onwards. Olf2_dur
        %will be listed as olf1_dur.
        
        dur_olf1_col_n = find_stim_mat_simple_col('duration', column_heads);        %identifying relevant column number in stim_mat_simple
        dur_olf2_col_n = find_stim_mat_simple_col('duration_olf2', column_heads);   %identifying relevant column number in stim_mat_simple
        
        raw_data_mat_orig = raw_data_mat;
        stim_mat_orig = stim_mat;
        
        %replacing NaNs with 0s in stim_mat_simple to do logic operations on
        %it's elements
        stim_mat_simple(isnan(stim_mat_simple)) = 0;
        
        %step1: Getting rid of all trials that aren't single-odor, single-pulse trials.
        db_od_trs = find(stim_mat_simple(:, dur_olf1_col_n) > 0.15 & stim_mat_simple(:, dur_olf2_col_n) > 0);   %dur_olf1 == 0.1 is just a place-holder, no stim actually delivered on olf_1.
        [raw_data_mat, stim_mat, stim_mat_simple, tr_list] = remove_trs(db_od_trs, raw_data_mat, stim_mat, stim_mat_simple, tr_list);         %removing trs with two odors delivered in the same trial
      
        
        %identifying multi-pulse train trials
        multi_pulse_trs = [];
        for tr_n = 1:size(stim_mat, 2)
            curr_train = stim_mat(tr_n).pulse_train;
            if size(curr_train, 1) > 1
                multi_pulse_trs = [multi_pulse_trs; tr_n];
            else
            end
        end
        [raw_data_mat, stim_mat, stim_mat_simple, tr_list] = remove_trs(multi_pulse_trs, raw_data_mat, stim_mat, stim_mat_simple, tr_list);   %removing trials with more than one long odor pulse
        
        
        %step2: Getting rid of all trials that aren't longest duration,
        %single pulse trials for all odors
        dur_list_olf1 = unique(stim_mat_simple(:, dur_olf1_col_n));
        dur_list_olf2 = unique(stim_mat_simple(:, dur_olf2_col_n));
        max_dur = max([dur_list_olf1; dur_list_olf2], [], 'omitnan');       %this is the longest duration stimulus, across olfactometers 
        short_dur_trs = find(stim_mat_simple(:, dur_olf1_col_n) < max_dur & stim_mat_simple(:, dur_olf2_col_n) < max_dur);
        [raw_data_mat, stim_mat, stim_mat_simple, tr_list] = remove_trs(short_dur_trs, raw_data_mat, stim_mat, stim_mat_simple, tr_list);   %removing trials with more shorter than max_dur odor pulses
        
        %combining odor numbers and durations for olf1 and olf2 so as not
        %to confuse the fitting program, which expects only one
        %olfactometer
        for tr_n = 1:size(stim_mat, 2)
            if isnan(stim_mat(tr_n).odours_olf2) == 0   %identifying an olf2 trial
                stim_mat(tr_n).odours = stim_mat(tr_n).odours_olf2 + 12;        %assigning an odor number that can later be decoded to be an olf2 odor if need be
                stim_mat(tr_n).duration = stim_mat(tr_n).duration_olf2;         %since multi-odor trials have been removed earlier, olf1_duration has to be 0.1 and can be discarded
            else
            end
        end
        %calculating dF/F traces from raw data
        filt_time = 0.2;            %in s, the time window for boxcar filter for generating filtered traces
        [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, curr_dir);
        [dff_data_mat_orig, dff_data_mat_f_orig] = cal_dff_traces_res(raw_data_mat_orig, stim_mat_orig, frame_time, filt_time, curr_dir);
        
        share_path = [share_path_base, dataset_list_name, '\fly', num2str(dir_n)];
        mkdir(share_path);
        save([share_path, '\dFF_data.mat'], 'dff_data_mat_f');
        save([share_path, '\stim_mat.mat'], 'stim_mat');
        save([share_path, '\path_orig.mat'], 'curr_dir');
        save([share_path, '\tr_list.mat'], 'tr_list');
        save([share_path, '\dff_data_mat_orig.mat'], 'dff_data_mat_f_orig');
        save([share_path, '\stim_mat_orig.mat'], 'stim_mat_orig');
        
    end
    
end

function [raw_data_mat, stim_mat, stim_mat_simple, tr_list] = remove_trs(rem_tr_list, raw_data_mat, stim_mat, stim_mat_simple, tr_list)
    raw_data_mat(:, :, rem_tr_list) = [];
    stim_mat_simple(rem_tr_list, :) = [];
    stim_mat(rem_tr_list) = [];
    tr_list(rem_tr_list) = [];
end