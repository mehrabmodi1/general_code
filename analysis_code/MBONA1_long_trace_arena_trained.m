clear all
close all

dir_list_path = 'C:\Data\Code\general_code_old\data_folder_lists\Janelia\Alpha1_60strace_71C03LxA_MB043CGal4_arena_trained.xls';

[del, dir_list] = xlsread(dir_list_path, 1);        %list of Suite2P results directories
n_dirs = size(dir_list, 1);
dataset_list_name = 'MBONA1_arena_paired_long_trace';

[del, odor_names1] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
[del, odor_names2] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList_olf2.xls', 1);

 %loop to go through all experiment datasets listed in list file
all_saved_resps = [];
fly_type_vec = [];
fly_n = 0;

skip_habit_trs = 1;     %1 - discards first two trials for each odor as odor habituation trials; 0 - includes all trials

for dir_n = 1:n_dirs
    fly_n = fly_n + 1;
    resp_mat = [];

    saved_an_results.scriptname = mfilename('fullpath');
    curr_dir = [dir_list{dir_n, 1}, '\'];
    curr_dir = manage_base_paths(curr_dir, 2);

    tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
    tif_times = tif_times.time_stamps;
    [stim_mat, stim_mat_simple, column_heads, color_vec, good_tr_list, params_orig] = load_params_trains_modular(curr_dir, tif_times);    %reading in trial stimulus parameters after matching time stamps to F tracesfore_colour = color_vec(1, :);
    distr_colour = color_vec(2, :);
    ctrl_colour = color_vec(3, :);

    %Reading in experimental parameters
    odor_list_olf1 = unique(stim_mat_simple(:, 1) );
    n_odors_olf1 = length(odor_list_olf1);
    odor_dur_list_olf1 = unique(stim_mat_simple(:, 2) );
    n_od_durs_olf1 = length(odor_dur_list_olf1);

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

    n_cells = size(raw_data_mat, 2);

    %calculating dF/F traces from raw data
    filt_time = 0.5;            %in s, the time window for boxcar filter for generating filtered traces
    [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, curr_dir);

    if size(dff_data_mat, 2) > 1
        dff_data_mat = mean(dff_data_mat, 2, 'omitnan');
        dff_data_mat_f = mean(dff_data_mat_f, 2, 'omitnan');
    else
    end
    %del = find(dff_data_mat_f < -1);
    %dff_data_mat_f(del) = -1;       %forcing crazy values to sane ones

    %identifying stim_mat_simple col numbers
    od_olf1_col_n = find_stim_mat_simple_col('odor_n', column_heads);           %identifying relevant column number in stim_mat_simple
    dur_olf1_col_n = find_stim_mat_simple_col('duration', column_heads);        %identifying relevant column number in stim_mat_simple
    od_durs = unique(stim_mat_simple(:, dur_olf1_col_n));
    od_durs(isnan(od_durs)) = [];
    
    %identifying fly type (left - 1 or right  - 0 arena) since paired odor would
    %differ by type
    if isempty(findstr(curr_dir, 'left')) == 0
        fly_type_vec = [fly_type_vec; 1];
    elseif isempty(findstr(curr_dir, 'right')) == 0
        fly_type_vec = [fly_type_vec; 0];
    else
        fly_type_vec = [fly_type_vec; nan];
    end
    
    saved_resp_mat = [];
    for od_ni = 1:length(odor_list_olf1)
        od_n = odor_list_olf1(od_ni);
        curr_trs = find(stim_mat_simple(:, od_olf1_col_n) == od_n);
        if skip_habit_trs == 1
            curr_trs(1:2) = [];
        else
        end
        
        
        curr_traces = squeeze(dff_data_mat_f(:, :, curr_trs));
        mean_trace = mean(curr_traces, 2, 'omitnan');
        
        %concatenating data across odors along dim 3
        saved_resp_mat = pad_n_concatenate(saved_resp_mat, curr_traces, 3, nan);
        
        figure(od_ni)
        %plotting resp traces
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
        stim_frs = stim_frs{1}; 
        plot(curr_traces, 'lineWidth', 2.5, 'Color', [0.65, 0.65, 0.65]);
        hold on
        plot(mean_trace, 'lineWidth', 2.5, 'Color', [0, 0, 0]);
        ylabel([odor_names1{od_n}, ' odor response (dF/F)']);
        set_xlabels_time(od_ni, frame_time, 10);
        fig_wrapup(od_ni, []);
        add_stim_bar(od_ni, stim_frs, color_vec(od_ni, :));
       
    end
    
    close figure 1
    close figure 2
    close figure 3
    
    %concatenating data across flies along dim 4
    all_saved_resps = pad_n_concatenate(all_saved_resps, saved_resp_mat, 4, nan);
    
    
    %PICK UP THREAD HERE
    %Keep track of fly types in fly_type_vec, quantify response sizes from all_saved_resps and
    %compare across fly types.
    
end