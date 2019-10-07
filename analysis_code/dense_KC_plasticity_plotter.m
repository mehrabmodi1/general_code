clear all
close all

dataset_list_paths = [...
                      
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_KC_dense_plasticity.xls'};...

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
        
        %In all flies, manually ensured that ROI 1 is the compartment with DAN innervation and 2, 3 have none.
        
        %plotting pre and post response traces for each compartment
        %paired odor, ROI 1
        figure(1)
        curr_trs = find(stim_mat_simple(:, od_olf1_col_n) == paired_odor_n);
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
        stim_frs = stim_frs{1};
        pre_trs = curr_trs(curr_trs < pairing_tr_n);
        post_trs = curr_trs(curr_trs > pairing_tr_n);
        pre_traces = squeeze(dff_data_mat_f(:, 1, pre_trs));
        post_traces = squeeze(dff_data_mat_f(:, 1, post_trs));        
        saved_responses(1, dir_n) = mean(mean(pre_traces(stim_frs(1):stim_frs(2) + round(4./frame_time), :), 1, 'omitnan'));
        saved_responses(2, dir_n) = mean(mean(post_traces(stim_frs(1):stim_frs(2) + round(4./frame_time), :), 1, 'omitnan'));
        
        plot(pre_traces, 'Color', color_vec(1, :), 'lineWidth', 2);
        hold on
        plot(post_traces, 'Color', color_vec(1, :).*0.6, 'lineWidth', 2);
        set_xlabels_time(1, frame_time, 25);
        ylabel('paired cpt, paired odor (dF/F)')
        fig_wrapup(1, script_name);
        add_stim_bar(1, stim_frs, [0.75, 0.75, 0.75]);
        
        
        %plotting pre and post response traces for each compartment
        %unpaired odor, ROI 1
        figure(2)
        curr_trs = find(stim_mat_simple(:, od_olf1_col_n) == unpaired_odor_n);
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
        stim_frs = stim_frs{1};
        pre_trs = curr_trs(curr_trs < pairing_tr_n);
        post_trs = curr_trs(curr_trs > pairing_tr_n);
        pre_traces = squeeze(dff_data_mat_f(:, 1, pre_trs));        
        post_traces = squeeze(dff_data_mat_f(:, 1, post_trs));        
        saved_responses(3, dir_n) = mean(mean(pre_traces(stim_frs(1):stim_frs(2) + round(4./frame_time), :), 1, 'omitnan'));
        saved_responses(4, dir_n) = mean(mean(post_traces(stim_frs(1):stim_frs(2) + round(4./frame_time), :), 1, 'omitnan'));
        
        plot(pre_traces, 'Color', color_vec(2, :), 'lineWidth', 2);
        hold on
        plot(post_traces, 'Color', color_vec(2, :).*0.6, 'lineWidth', 2);
        set_xlabels_time(2, frame_time, 25);
        ylabel('paired cpt, unpaired odor (dF/F)')
        fig_wrapup(2, script_name);
        add_stim_bar(2, stim_frs, [0.75, 0.75, 0.75]);
        
        
        %plotting pre and post response traces for each compartment
        %paired odor, ROI 2
        figure(3)
        curr_trs = find(stim_mat_simple(:, od_olf1_col_n) == paired_odor_n);
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
        stim_frs = stim_frs{1};
        pre_trs = curr_trs(curr_trs < pairing_tr_n);
        post_trs = curr_trs(curr_trs > pairing_tr_n);
        pre_traces = squeeze(dff_data_mat_f(:, 2, pre_trs));        
        post_traces = squeeze(dff_data_mat_f(:, 2, post_trs));
        saved_responses(5, dir_n) = mean(mean(pre_traces(stim_frs(1):stim_frs(2) + round(4./frame_time), :), 1, 'omitnan'));
        saved_responses(6, dir_n) = mean(mean(post_traces(stim_frs(1):stim_frs(2) + round(4./frame_time), :), 1, 'omitnan'));
        
        plot(pre_traces, 'Color', color_vec(1, :), 'lineWidth', 2);
        hold on
        plot(post_traces, 'Color', color_vec(1, :).*0.6, 'lineWidth', 2);
        set_xlabels_time(3, frame_time, 25);
        ylabel('unpaired cpt, paired odor (dF/F)')
        fig_wrapup(3, script_name);
        add_stim_bar(3, stim_frs, [0.75, 0.75, 0.75]);
        
        
        %plotting pre and post response traces for each compartment
        %unpaired odor, ROI 2
        figure(4)
        curr_trs = find(stim_mat_simple(:, od_olf1_col_n) == unpaired_odor_n);
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
        stim_frs = stim_frs{1};
        pre_trs = curr_trs(curr_trs < pairing_tr_n);
        post_trs = curr_trs(curr_trs > pairing_tr_n);
        pre_traces = squeeze(dff_data_mat_f(:, 1, pre_trs));        
        post_traces = squeeze(dff_data_mat_f(:, 1, post_trs));
        saved_responses(7, dir_n) = mean(mean(pre_traces(stim_frs(1):stim_frs(2) + round(4./frame_time), :), 1, 'omitnan'));
        saved_responses(8, dir_n) = mean(mean(post_traces(stim_frs(1):stim_frs(2) + round(4./frame_time), :), 1, 'omitnan'));
        
        plot(pre_traces, 'Color', color_vec(2, :), 'lineWidth', 2);
        hold on
        plot(post_traces, 'Color', color_vec(2, :).*0.6, 'lineWidth', 2);
        set_xlabels_time(4, frame_time, 25);
        ylabel('unpaired cpt, unpaired odor (dF/F)')
        fig_wrapup(4, script_name);
        add_stim_bar(4, stim_frs, [0.75, 0.75, 0.75]);
        
        
        keyboard
        
        close figure 1
        close figure 2
        close figure 3
        close figure 4
                
    end
    marker_colors = [color_vec(1, :); color_vec(1, :).*0.6; color_vec(1, :); color_vec(1, :).*0.6; color_vec(2, :); color_vec(2, :).*0.6; color_vec(2, :); color_vec(2, :).*0.6 ];
    col_pairs = [1, 2; 3, 4; 5, 6; 7, 8];
    scattered_dot_plot(saved_responses', 5, 1, 4, 8, marker_colors, 1, col_pairs, [0.75, 0.75, 0.75],...
                            [{'yellow-paired_p_r_e'}, {'yellow-paired_p_s_t'}, {'yellow-unpaired_p_r_e'}, {'yellow-unpaired_p_s_t'},...
                                {'green-paired_p_r_e'}, {'green-paired_p_s_t'}, {'green-unpaired_p_r_e'}, {'green-unpaired_p_s_t'}], 1, [0.35, 0.35, 0.35]);

end