clear all
close all

dataset_list_paths = [...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_Alpha1_60strace_71C03LxA_MB043CGal4.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_Alpha1_60strace_72D01LxAChr88tdT.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_Alpha1_60strace_71C03LxA_MB043CGal4_noChrisoncontrol.xls'};...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\Alpha1_60strace_71C03LxA_MB043CGal4_Chrison_noLED_control.xls'};...                      
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
   
    %loop to go through all experiment datasets listed in list file
    saved_resps = [];
    for dir_n = 1:n_dirs
        fly_n = fly_n + 1;
              
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
        led_on_col_n = find_stim_mat_simple_col('led_on', column_heads);            %identifying relevant column number in stim_mat_simple
        od_olf1_col_n = find_stim_mat_simple_col('odor_n', column_heads);           %identifying relevant column number in stim_mat_simple
        dur_olf1_col_n = find_stim_mat_simple_col('duration', column_heads);        %identifying relevant column number in stim_mat_simple
        od_durs = unique(stim_mat_simple(:, dur_olf1_col_n));
        od_durs(isnan(od_durs)) = [];
                
        %identifying relevant odor numbers for each olfactometer
        pairing_trs = find(stim_mat_simple(:, led_on_col_n) == 1);
        if isempty(pairing_trs) == 1
            pairing_trs = (((size(dff_data_mat, 3) - 6)./2) + 1):((size(dff_data_mat, 3) - 6)./2) + 7;
        else
        end
        paired_od_n = stim_mat_simple(pairing_trs(1), od_olf1_col_n);
        unpaired_od_n = [3, 11];        %only PA or EL are ever paired with LED
        unpaired_od_n(unpaired_od_n == paired_od_n) = [];
        ctrl_od_n = 9;
        
        pre_trs = 1:1:(min(pairing_trs) - 1);
        post_trs = (max(pairing_trs) + 2):1:size(dff_data_mat, 3);          %accounting for the CS- trial after the last pairing trial
        
        
        %plotting responses for the paired odor
        curr_trs_pre = find(stim_mat_simple(:, od_olf1_col_n) == paired_od_n);
        curr_trs_pre = intersect(curr_trs_pre, pre_trs);
        %curr_trs_pre = curr_trs_pre(1:7);                             %discarding first presentation of any odor for habituation effects
        curr_trs_post = find(stim_mat_simple(:, od_olf1_col_n) == paired_od_n);
        curr_trs_post = intersect(curr_trs_post, post_trs);
        %curr_trs_post = curr_trs_post(1:7);  
        pre_traces = squeeze(dff_data_mat_f(:, 1, curr_trs_pre));
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs_pre(1), frame_time);
        stim_frs = stim_frs{1};      
        
        mean_pre_trace = mean(pre_traces, 2);
        resp_vec(1, 1) = max(mean_pre_trace(stim_frs(1):(stim_frs(2) + round(2./frame_time))));       %taking max of filtered trace as resp measure
        post_traces = squeeze(dff_data_mat_f(:, 1, curr_trs_post));
        mean_post_trace = mean(post_traces, 2);
        resp_vec(1, 2) = max(mean_post_trace(stim_frs(1):(stim_frs(2) + round(2./frame_time))));      %taking max of filtered trace as resp measure
                
                    
        
        figure(1)
        plot(pre_traces, 'Color', [0.6, 0.6, 0.6], 'lineWidth', 1.5)
        hold on
        plot(post_traces, 'Color', [0.1, 0.1, 0.1], 'lineWidth', 1.5)
        ylabel([odor_names1{paired_od_n}, ' odor response (dF/F)']);
        set_xlabels_time(1, frame_time, 10);
        fig_wrapup(1, []);
        add_stim_bar(1, stim_frs, color_vec(2, :));
        
        
        
        %plotting responses for the un-paired odor
        curr_trs_pre = find(stim_mat_simple(:, od_olf1_col_n) == unpaired_od_n);
        curr_trs_pre = intersect(curr_trs_pre, pre_trs);
        %curr_trs_pre = curr_trs_pre(1:7);                           %discarding first presentation of any odor for habituation effects
        curr_trs_post = find(stim_mat_simple(:, od_olf1_col_n) == unpaired_od_n);
        curr_trs_post = intersect(curr_trs_post, post_trs);
        %curr_trs_post = curr_trs_post(1:7);
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs_pre(1), frame_time);
        stim_frs = stim_frs{1};      
        
        pre_traces = squeeze(dff_data_mat_f(:, 1, curr_trs_pre));
        mean_pre_trace = mean(pre_traces, 2);
        resp_vec(1, 3) = max(mean_pre_trace(stim_frs(1):(stim_frs(2) + round(2./frame_time))));       %taking max of filtered trace as resp measure
        post_traces = squeeze(dff_data_mat_f(:, 1, curr_trs_post));
        mean_post_trace = mean(post_traces, 2);
        resp_vec(1, 4) = max(mean_post_trace(stim_frs(1):(stim_frs(2) + round(2./frame_time))));      %taking max of filtered trace as resp measure
        
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs_pre(1), frame_time);
        stim_frs = stim_frs{1};                   
        
        figure(2)
        plot(pre_traces, 'Color', [0.6, 0.6, 0.6], 'lineWidth', 1.5)
        hold on
        plot(post_traces, 'Color', [0.1, 0.1, 0.1], 'lineWidth', 1.5)
        ylabel([odor_names1{unpaired_od_n}, ' odor response (dF/F)']);
        set_xlabels_time(2, frame_time, 10);
        fig_wrapup(2, []);
        add_stim_bar(2, stim_frs, color_vec(1, :));
        
        
        %plotting responses for the third, ctrl odor
        curr_trs_pre = find(stim_mat_simple(:, od_olf1_col_n) == ctrl_od_n);
        curr_trs_pre = intersect(curr_trs_pre, pre_trs);
        %curr_trs_pre = curr_trs_pre(1:7);                           %discarding first presentation of any odor for habituation effects
        curr_trs_post = find(stim_mat_simple(:, od_olf1_col_n) == ctrl_od_n);
        curr_trs_post = intersect(curr_trs_post, post_trs);
        %curr_trs_post = curr_trs_post(1:7);
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs_pre(1), frame_time);
        stim_frs = stim_frs{1};      
        
        pre_traces = squeeze(dff_data_mat_f(:, 1, curr_trs_pre));
        mean_pre_trace = mean(pre_traces, 2);
        resp_vec(1, 3) = max(mean_pre_trace(stim_frs(1):(stim_frs(2) + round(2./frame_time))));       %taking max of filtered trace as resp measure
        post_traces = squeeze(dff_data_mat_f(:, 1, curr_trs_post));
        mean_post_trace = mean(post_traces, 2);
        resp_vec(1, 4) = max(mean_post_trace(stim_frs(1):(stim_frs(2) + round(2./frame_time))));      %taking max of filtered trace as resp measure
        
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs_pre(1), frame_time);
        stim_frs = stim_frs{1};                   
        
        figure(3)
        plot(pre_traces, 'Color', [0.6, 0.6, 0.6], 'lineWidth', 1.5)
        hold on
        plot(post_traces, 'Color', [0.1, 0.1, 0.1], 'lineWidth', 1.5)
        ylabel('control odor response (dF/F)');
        set_xlabels_time(3, frame_time, 10);
        fig_wrapup(3, []);
        add_stim_bar(3, stim_frs, color_vec(3, :));
        
        
        
        if suppress_plots == 0
            keyboard
        else
        end
        
        close figure 1
        close figure 2
        close figure 3
        
        saved_resps = [saved_resps; resp_vec];
       
    end
    marker_colors = [color_vec(1, :); color_vec(1, :); color_vec(2, :); color_vec(2, :)];
    col_pairs = [1, 2; 3, 4];
    scattered_dot_plot(saved_resps(:, 1:4), 3, 1, 4, 8, marker_colors, 1, col_pairs, [0.75, 0.75, 0.75],...
                            [{'paired_p_r_e'}, {'paired_p_o_s_t'}, {'unpaired_p_r_e'}, {'unpaired_p_o_s_t'}], 1, [0.35, 0.35, 0.35]);
    
    ylabel('response size (mean dF/F)')
    fig_wrapup(3, []);
    
    %statistical testing
    [h_paired, p_paired] = ttest(saved_resps(:, 1), saved_resps(:, 2))
    [h_unpaired, p_unpaired] = ttest(saved_resps(:, 3), saved_resps(:, 4))
    [h_post_comp, p_post_comp] = ttest(saved_resps(:, 2), saved_resps(:, 4))
    [h_ratios, p_ratios] = ttest(saved_resps(:, 1)./saved_resps(:, 2), saved_resps(:, 3)./saved_resps(:, 4))

end