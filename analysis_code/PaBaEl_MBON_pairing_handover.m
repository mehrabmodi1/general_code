clear all
close all

dataset_list_paths = [...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONAlpha1.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONAlpha1_set2.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONAlpha1_set3_highLED.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONGamma2.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONGamma2_set1.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONGamma2_set2.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_PaBaEl_MBONG2_handover.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_MBONG2_PaBaEl_handover_simple_starved.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_MBONG2_PaBaEl_handover_starved36_halfAra.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_MBONG2_PaBaEl_handover_starved_halfAra.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_MBONG2_PaBaEl_handover_starved_halfAra_prehabituated.xls'};...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_MBONG2_PaBaEl_handover_starved_halfAra_prehabituated_strongUS.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_MBONG2_PaBaEl_handover_starved_halfAra_prehabituated_strongUS_EL_handover.xls'};...
                      
                      ];
            
suppress_plots = 0;
[del, odor_names1] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
[del, odor_names2] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList_olf2.xls', 1);
odor_names2{3} = 'Butyl acetate';

global color_vec;                
a = colormap('bone');
global greymap
greymap = flipud(a);
fly_n = 0;
script_name = mfilename;

y_ax_traces = 0.8;
y_ax_fit_traces = 0.6;

for list_n = 1:size(dataset_list_paths, 1)
    curr_dir_list_path = dataset_list_paths{list_n, 1};
    [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
    n_dirs = size(dir_list, 1);
    dataset_list_name = findstr(curr_dir_list_path, 'list_');
    dataset_list_name = curr_dir_list_path((dataset_list_name + 5):(end - 4));
    
    
    %loop to go through all experiment datasets listed in list file
    saved_resps_hover = [];
    saved_resps_simple = [];
    saved_mean_traces_simple = [];
    saved_mean_traces_transition = [];
    
    
    for dir_n = 1:n_dirs
        fly_n = fly_n + 1;
              
        saved_an_results.scriptname = mfilename('fullpath');
        curr_dir = [dir_list{dir_n, 1}, '\'];
        curr_dir = manage_base_paths(curr_dir, 2);
       
        tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
        tif_times = tif_times.time_stamps;
        [stim_mat, stim_mat_simple, column_heads, color_vec, good_tr_list, params_orig] = load_params_trains_modular(curr_dir, tif_times);    %reading in trial stimulus parameters after matching time stamps to F traces
        paired_color = color_vec(2, :);
        unpaired_color = color_vec(1, :);
        EL_color = color_vec(3, :);
        
        %Reading in experimental parameters
        odor_list_olf1 = unique(stim_mat_simple(:, 1) );
        n_odors_olf1 = length(odor_list_olf1);
        odor_dur_list_olf1 = unique(stim_mat_simple(:, 2) );
        n_od_durs_olf1 = length(odor_dur_list_olf1);
        n_trains_olf1 = max(stim_mat_simple(:, 8));
        
        odor_list_olf2 = unique(stim_mat_simple(:, 9) );
        n_odors_olf2 = length(odor_list_olf2);
        odor_dur_list_olf2 = unique(stim_mat_simple(:, 10) );
        n_od_durs_olf2 = length(odor_dur_list_olf2);
        n_trains_olf2 = max(stim_mat_simple(:, 17));
        
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
        
        %dumping data from manually identified, z-drifted trials
        bad_tr_list = 1:size(raw_data_mat, 3);
        bad_tr_list(good_tr_list) = [];
        raw_data_mat(:, :, bad_tr_list) = nan;
        
        n_cells = size(raw_data_mat, 2);
        
        %calculating dF/F traces from raw data
        filt_time = 0.2;            %in s, the time window for boxcar filter for generating filtered traces
        [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, curr_dir);
        
        if size(dff_data_mat, 2) > 1
            dff_data_mat = mean(dff_data_mat, 2, 'omitnan');
            dff_data_mat_f = mean(dff_data_mat_f, 2, 'omitnan');
        else
        end
        del = find(dff_data_mat_f < -1);
        dff_data_mat_f(del) = -1;       %forcing crazy values to sane ones
        
        %identifying stim_mat_simple col numbers
        led_on_col_n = find_stim_mat_simple_col('led_on', column_heads);            %identifying relevant column number in stim_mat_simple
        od_olf1_col_n = find_stim_mat_simple_col('odor_n', column_heads);           %identifying relevant column number in stim_mat_simple
        od_olf2_col_n = find_stim_mat_simple_col('odour_olf2', column_heads);       %identifying relevant column number in stim_mat_simple
        dur_olf1_col_n = find_stim_mat_simple_col('duration', column_heads);        %identifying relevant column number in stim_mat_simple
        dur_olf2_col_n = find_stim_mat_simple_col('duration_olf2', column_heads);   %identifying relevant column number in stim_mat_simple
        od_col_ns = [od_olf1_col_n, od_olf2_col_n];
        dur_col_ns = [dur_olf1_col_n, dur_olf2_col_n];
        od_durs = unique(stim_mat_simple(:, dur_col_ns(2)));
        od_durs(isnan(od_durs)) = [];
        odn_list_olf2 = unique(stim_mat_simple(:, od_col_ns(2)));
        odn_list_olf2(isnan(odn_list_olf2)) = [];
        
        
        %identifying relevant odor numbers for each olfactometer
        pairing_tr_n = find(stim_mat_simple(:, led_on_col_n) == 1);        
        
        paired_od_n_olf2 = stim_mat_simple(pairing_tr_n, od_olf2_col_n);    %paired odor is always an olf2 odor        
        paired_od_name = odor_names2{paired_od_n_olf2};
        paired_od_n_olf1 = od_name_lookup(odor_names1, paired_od_name); 

        %unpaired_od_n_olf2 = odn_list_olf2;
        %unpaired_od_n_olf2(unpaired_od_n_olf2 == intersect(odn_list_olf2, paired_od_n_olf2)) = [];

        %if paired odor is PA, unpaired odor must be BA or vice versa
        unpaired_od_n_olf2 = [1, 3];
        [del, deli] = intersect(unpaired_od_n_olf2, paired_od_n_olf2);
        unpaired_od_n_olf2(deli)= [];

        unpaired_od_name = odor_names2{unpaired_od_n_olf2};
        unpaired_od_n_olf1 = od_name_lookup(odor_names1, unpaired_od_name); 
        
        %list of all handover trials, with any combination of odors
        handover_trs = find(stim_mat_simple(:, dur_col_ns(1)) == od_durs(1) & stim_mat_simple(:, dur_col_ns(2)) == od_durs(1));
        
        %building tr lists
        all_trs_paired_od_olf2 = find(stim_mat_simple(:, od_col_ns(2)) == paired_od_n_olf2 & stim_mat_simple(:, dur_col_ns(2)) == od_durs(1));   %this includes post-pairing and handover trials
        all_trs_unpaired_od_olf1 = find(stim_mat_simple(:, od_col_ns(1)) == unpaired_od_n_olf1 & stim_mat_simple(:, dur_col_ns(1)) == od_durs(1));   %this includes post-pairing and handover trials
        all_trs_paired_od_olf1 = find(stim_mat_simple(:, od_col_ns(1)) == paired_od_n_olf1 & stim_mat_simple(:, dur_col_ns(1)) == od_durs(1));   %this includes post-pairing and handover trials
        all_trs_unpaired_od_olf2 = find(stim_mat_simple(:, od_col_ns(2)) == unpaired_od_n_olf2 & stim_mat_simple(:, dur_col_ns(2)) == od_durs(1));   %this includes post-pairing and handover trials
        
        
        %1. non-handover, paired odor
        curr_trs = union(all_trs_paired_od_olf2, all_trs_paired_od_olf1);
        [del, del2] = intersect(curr_trs, handover_trs);
        curr_trs(del2) = [];            %excluded handover trs
        curr_trs_pre = curr_trs(curr_trs < pairing_tr_n);   %only pre trs 
        curr_trs_post = curr_trs(curr_trs > pairing_tr_n);   %only post trs 
                
        curr_traces_pre = dff_data_mat_f(:, :, curr_trs_pre(2:end));
        mean_trace_pre = mean(squeeze(curr_traces_pre), 2, 'omitnan');
        curr_traces_post = dff_data_mat_f(:, :, curr_trs_post(2:end));
        mean_trace_post = mean(squeeze(curr_traces_post), 2, 'omitnan');
        paired_pre_mean = mean(squeeze(curr_traces_pre), 2, 'omitnan');
        paired_post_mean = mean(squeeze(curr_traces_post), 2, 'omitnan');
        length_data_paired = length(paired_post_mean) - sum(isnan(paired_post_mean));
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs_pre(1), frame_time);
        stim_frs = stim_frs{2};                   %because olf2 is used for paired odor
        
        paired_simple_pre = max(mean(squeeze(curr_traces_pre(stim_frs(1, 1):(stim_frs(1, 2)+ round(3./frame_time) ), :)), 2, 'omitnan'));
        paired_simple_post = max(mean(squeeze(curr_traces_post(stim_frs(1, 1):(stim_frs(1, 2)+ round(3./frame_time) ), :)), 2, 'omitnan'));
        
        paired_simple_trace_pre = mean(squeeze(curr_traces_pre(stim_frs(1, 1):(stim_frs(1, 2)+ round(3./frame_time) ), :)), 2, 'omitnan');
        paired_simple_trace_post = mean(squeeze(curr_traces_post(stim_frs(1, 1):(stim_frs(1, 2)+ round(3./frame_time) ), :)), 2, 'omitnan');
        
        
        figure(1)
        plot(squeeze(curr_traces_pre), 'Color', [0.6, 0.6, 0.6], 'lineWidth', 0.5);
        hold on
        plot(squeeze(curr_traces_post), 'Color', [0.1, 0.1, 0.1], 'lineWidth', 0.5);
        plot(mean_trace_pre, 'Color', [0.6, 0.6, 0.6], 'lineWidth', 2);
        plot(mean_trace_post, 'Color', [0.1, 0.1, 0.1], 'lineWidth', 2);
        
        %shadedErrorBar([], trace_mean, trace_se, {'Color', fore_colour})
        ylabel('paired odor, simple (dF/F)')
        ax_vals = axis;
        ax_vals(1, 4) = y_ax_traces;
        axis(ax_vals);
        set_xlabels_time(1, frame_time, 10);
        fig_wrapup(1, script_name);
        add_stim_bar(1, stim_frs, paired_color);
        hold off
        
        
        %3. non-handover, unpaired odor
        curr_trs = union(all_trs_unpaired_od_olf1, all_trs_unpaired_od_olf2);
        [del, del2] = intersect(curr_trs, handover_trs);
        curr_trs(del2) = [];            %excluded handover trs
        curr_trs_pre = curr_trs(curr_trs < pairing_tr_n);   %only pre trs 
        curr_trs_post = curr_trs(curr_trs > pairing_tr_n);   %only post trs 
                
        curr_traces_pre = dff_data_mat_f(:, :, curr_trs_pre(2:end));
        mean_trace_pre = mean(squeeze(curr_traces_pre), 2, 'omitnan');
        curr_traces_post = dff_data_mat_f(:, :, curr_trs_post(2:end));
        mean_trace_post = mean(squeeze(curr_traces_post), 2, 'omitnan');
        
        unpaired_pre_mean = mean(squeeze(curr_traces_pre), 2);
        unpaired_post_mean = mean(squeeze(curr_traces_pre), 2);
        length_data_unpaired = length(unpaired_post_mean) - sum(isnan(unpaired_post_mean));
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs_pre(1), frame_time);
        stim_frs = stim_frs{1};                   %because olf2 is used for paired odor
        
        unpaired_simple_pre = max(mean(squeeze(curr_traces_pre(stim_frs(1, 1):(stim_frs(1, 2)+ round(3./frame_time) ), :)), 2, 'omitnan'));
        unpaired_simple_post = max(mean(squeeze(curr_traces_post(stim_frs(1, 1):(stim_frs(1, 2)+ round(3./frame_time) ), :)), 2, 'omitnan'));
        
        unpaired_simple_trace_pre = mean(squeeze(curr_traces_pre(stim_frs(1, 1):(stim_frs(1, 2)+ round(3./frame_time) ), :)), 2, 'omitnan');
        unpaired_simple_trace_post = mean(squeeze(curr_traces_post(stim_frs(1, 1):(stim_frs(1, 2)+ round(3./frame_time) ), :)), 2, 'omitnan');
        
        
        figure(3)
        plot(squeeze(curr_traces_pre), 'Color', [0.6, 0.6, 0.6], 'lineWidth', 0.5);
        hold on
        plot(squeeze(curr_traces_post), 'Color', [0.1, 0.1, 0.1], 'lineWidth', 0.5);
        plot(mean_trace_pre, 'Color', [0.6, 0.6, 0.6], 'lineWidth', 2);
        plot(mean_trace_post, 'Color', [0.1, 0.1, 0.1], 'lineWidth', 2);
        ylabel('unpaired odor, simple (dF/F)')
        ax_vals = axis;
        ax_vals(1, 4) = y_ax_traces;
        axis(ax_vals);
        set_xlabels_time(3, frame_time, 10);
        fig_wrapup(3, script_name);
        add_stim_bar(3, stim_frs, unpaired_color);
        hold off
        
        %getting rid of nans
        length_real_data = min([length_data_paired; length_data_unpaired]);
        paired_pre_mean( (length_real_data + 1):end) = [];
        paired_post_mean( (length_real_data + 1):end) = [];
        unpaired_pre_mean( (length_real_data + 1):end) = [];
        unpaired_post_mean( (length_real_data + 1):end) = [];
        
        %computing linear sums with time-offset to simulate handover trials
        n_fr_delay = round(od_durs(1)./frame_time);
        pad = zeros(n_fr_delay, 1) + nan;
        plusminus_sum_pre = [paired_pre_mean; pad] + [pad; unpaired_pre_mean];           %CS+ first, CS- second
        plusminus_sum_post = [paired_post_mean; pad] + [pad; unpaired_post_mean];        %CS+ first, CS- second
        
        minusplus_sum_pre = sum([[pad; paired_pre_mean], [unpaired_pre_mean; pad]], 2);          %CS- first, CS + second
        minusplus_sum_post = sum([[pad; paired_post_mean], [unpaired_post_mean; pad]], 2);          %CS- first, CS + second

        
        %2. handover, paired odor
        curr_trs = handover_trs;
        curr_trs = intersect(curr_trs, all_trs_paired_od_olf2);
        if stim_mat_simple(curr_trs(1), od_col_ns(1)) == 11      %ie. EL on olf1
            handover_color1 = EL_color;
        else
            handover_color1 = unpaired_color;
        end
        curr_trs_pre = curr_trs(curr_trs < pairing_tr_n);   %only pre trs 
        curr_trs_post = curr_trs(curr_trs > pairing_tr_n);   %only post trs 
                
        curr_traces_pre = dff_data_mat_f(:, :, curr_trs_pre(2:end));
        mean_trace_pre = mean(squeeze(curr_traces_pre), 2, 'omitnan');
        curr_traces_post = dff_data_mat_f(:, :, curr_trs_post(2:end));
        mean_trace_post = mean(squeeze(curr_traces_post), 2, 'omitnan');
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs_pre(1), frame_time);
        stim_frs = [stim_frs{1}; stim_frs{2}];                   %because olf2 is used for paired odor
        
        %computing response to second odor pulse in handover presentations
        paired_mean_resp_pre = max(mean(squeeze(curr_traces_pre(stim_frs(2, 1):(stim_frs(2, 2)+ round(3./frame_time) ), :)), 2, 'omitnan'));
        paired_mean_resp_post = max(mean(squeeze(curr_traces_post(stim_frs(2, 1):(stim_frs(2, 2)+ round(3./frame_time) ), :)), 2, 'omitnan'));
        
        paired_trans_trace_pre = mean(squeeze(curr_traces_pre(stim_frs(1, 1):(stim_frs(1, 2)+ round(3./frame_time) ), :)), 2, 'omitnan');
        paired_trans_trace_post = mean(squeeze(curr_traces_post(stim_frs(1, 1):(stim_frs(1, 2)+ round(3./frame_time) ), :)), 2, 'omitnan');
        
        
        %raw trace plot
        figure(2)
        plot(squeeze(curr_traces_pre), 'Color', [0.6, 0.6, 0.6], 'lineWidth', 0.5);
        hold on
        plot(squeeze(curr_traces_post), 'Color', [0.1, 0.1, 0.1], 'lineWidth', 0.5);
        plot(mean_trace_pre, 'Color', [0.6, 0.6, 0.6], 'lineWidth', 2);
        plot(mean_trace_post, 'Color', [0.1, 0.1, 0.1], 'lineWidth', 2);
        ylabel('paired odor, handover (dF/F)')
        ax_vals = axis;
        ax_vals(1, 4) = y_ax_traces;
        axis(ax_vals);
        set_xlabels_time(2, frame_time, 10);
        fig_wrapup(2, script_name);
        add_stim_bar(2, stim_frs, [handover_color1; paired_color]);
        hold off
        
        
        %linear model plot
        figure(5)
        mean_pre = mean(squeeze(curr_traces_pre), 2, 'omitnan');
        se_pre = std(squeeze(curr_traces_pre), [], 2, 'omitnan');%./sqrt(size(squeeze(curr_traces_pre), 2));
        mean_post = mean(squeeze(curr_traces_post), 2, 'omitnan');
        se_post = std(squeeze(curr_traces_post), [], 2, 'omitnan');%./sqrt(size(squeeze(curr_traces_post), 2));
        shadedErrorBar([], mean_pre, se_pre, {'Color', [0.6, 0.6, 0.6], 'lineWidth', 1.5});
        hold on
        shadedErrorBar([], mean_post, se_post, {'Color', [0, 0, 0], 'lineWidth', 1.5});
        plot(minusplus_sum_pre, 'Color', [0.5, 0.7, 0.9], 'lineWidth', 1.5);
        plot(minusplus_sum_post, 'Color', [0.2, 0.35, 0.6], 'lineWidth', 1.5);
        ylabel('paired odor, handover (dF/F)')
        ax_vals = axis;
        ax_vals(1, 4) = y_ax_fit_traces;
        axis(ax_vals);
        set_xlabels_time(5, frame_time, 10);
        fig_wrapup(5, script_name);
        add_stim_bar(5, stim_frs, [handover_color1; paired_color;]);
        hold off
        
        %4. handover, unpaired odor
        curr_trs = handover_trs;
        curr_trs = intersect(curr_trs, all_trs_unpaired_od_olf2);
        curr_trs_pre = curr_trs(curr_trs < pairing_tr_n);   %only pre trs 
        curr_trs_post = curr_trs(curr_trs > pairing_tr_n);   %only post trs 
                
        if stim_mat_simple(curr_trs(1), od_col_ns(1)) == 11      %ie. EL on olf1
            handover_color1 = EL_color;
        elseif stim_mat_simple(curr_trs(1), od_col_ns(1)) == paired_od_n_olf1
            handover_color1 = paired_color;
        elseif stim_mat_simple(curr_trs(1), od_col_ns(1)) == unpaired_od_n_olf1
            handover_color1 = unpaired_color;
        end
        
        curr_traces_pre = dff_data_mat_f(:, :, curr_trs_pre(2:end));
        curr_traces_post = dff_data_mat_f(:, :, curr_trs_post(2:end));
        mean_trace_pre = mean(squeeze(curr_traces_pre), 2, 'omitnan');
        mean_trace_post = mean(squeeze(curr_traces_post), 2, 'omitnan');
        
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs_pre(1), frame_time);
        stim_frs = [stim_frs{1}; stim_frs{2}];                   %because olf2 is used for paired odor
        
        %computing response to second odor pulse in handover presentations
        unpaired_mean_resp_pre = max(mean(squeeze(curr_traces_pre(stim_frs(2, 1):(stim_frs(2, 2) + round(3./frame_time) ), :)), 2, 'omitnan'));
        unpaired_mean_resp_post = max(mean(squeeze(curr_traces_post(stim_frs(2, 1):(stim_frs(2, 2)+ round(3./frame_time) ), :)), 2, 'omitnan'));
        
        unpaired_trans_trace_pre = mean(squeeze(curr_traces_pre(stim_frs(1, 1):(stim_frs(1, 2)+ round(3./frame_time) ), :)), 2, 'omitnan');
        unpaired_trans_trace_post = mean(squeeze(curr_traces_post(stim_frs(1, 1):(stim_frs(1, 2)+ round(3./frame_time) ), :)), 2, 'omitnan');
        
        
        saved_resps_hover = [saved_resps_hover; paired_mean_resp_pre, paired_mean_resp_post, unpaired_mean_resp_pre, unpaired_mean_resp_post];
        saved_resps_simple = [saved_resps_simple; paired_simple_pre, paired_simple_post, unpaired_simple_pre, unpaired_simple_post];
        
        saved_mean_traces_simple = cat(3, saved_mean_traces_simple, [paired_simple_trace_pre, paired_simple_trace_post, unpaired_simple_trace_pre, unpaired_simple_trace_post]);
        saved_mean_traces_transition = cat(3, saved_mean_traces_transition, [paired_trans_trace_pre, paired_trans_trace_post, unpaired_trans_trace_pre, unpaired_trans_trace_post]);
        
        figure(4)
        plot(squeeze(curr_traces_pre), 'Color', [0.6, 0.6, 0.6], 'lineWidth', 0.5);
        hold on
        plot(squeeze(curr_traces_post), 'Color', [0.1, 0.1, 0.1], 'lineWidth', 0.5);
        plot(mean_trace_pre, 'Color', [0.6, 0.6, 0.6], 'lineWidth', 2);
        plot(mean_trace_post, 'Color', [0.1, 0.1, 0.1], 'lineWidth', 2);
        ylabel('unpaired odor, handover (dF/F)')
        ax_vals = axis;
        ax_vals(1, 4) = y_ax_traces;
        axis(ax_vals);
        set_xlabels_time(4, frame_time, 10);
        fig_wrapup(4, script_name);
        add_stim_bar(4, stim_frs, [handover_color1; unpaired_color;]);
        hold off
        
        
        %linear model plot
        figure(6)
        mean_pre = mean(squeeze(curr_traces_pre), 2, 'omitnan');
        se_pre = std(squeeze(curr_traces_pre), [], 2, 'omitnan');%./sqrt(size(squeeze(curr_traces_pre), 2));
        mean_post = mean(squeeze(curr_traces_post), 2, 'omitnan');
        se_post = std(squeeze(curr_traces_post), [], 2, 'omitnan');%./sqrt(size(squeeze(curr_traces_post), 2));
        shadedErrorBar([], mean_pre, se_pre, {'Color', [0.6, 0.6, 0.6], 'lineWidth', 1.5});
        hold on
        shadedErrorBar([], mean_post, se_post, {'Color', [0, 0, 0], 'lineWidth', 1.5});
        plot(plusminus_sum_pre, 'Color', [0.5, 0.7, 0.9], 'lineWidth', 1.5);
        plot(plusminus_sum_post, 'Color', [0.2, 0.35, 0.6], 'lineWidth', 1.5);
        ylabel('unpaired odor, handover (dF/F)')
        ax_vals = axis;
        ax_vals(1, 4) = y_ax_fit_traces;
        axis(ax_vals);
        set_xlabels_time(6, frame_time, 10);
        fig_wrapup(6, script_name);
        add_stim_bar(6, stim_frs, [handover_color1; unpaired_color;]);
        hold off
        
        
        
        if suppress_plots == 0
            keyboard
        else
        end
      
     
        close figure 1
        close figure 2
        close figure 3
        close figure 4
        close figure 5
        close figure 6
                
    end
    marker_colors = [paired_color; paired_color; unpaired_color; unpaired_color];
    col_pairs = [1, 2; 3, 4];
    scattered_dot_plot(saved_resps_hover(:, 1:4), 7, 1, 4, 8, marker_colors, 1, col_pairs, [0.75, 0.75, 0.75],...
                            [{'paired_p_r_e'}, {'paired_p_o_s_t'}, {'unpaired_p_r_e'}, {'unpaired_p_o_s_t'}], 1, [0.35, 0.35, 0.35]);
    
    ylabel('second pulse resp (mean dF/F)')
    fig_wrapup(7, script_name);
    
    
    
    marker_colors = [paired_color; paired_color; unpaired_color; unpaired_color];
    col_pairs = [1, 2; 3, 4];
    scattered_dot_plot(saved_resps_simple(:, 1:4), 8, 1, 4, 8, marker_colors, 1, col_pairs, [0.75, 0.75, 0.75],...
                            [{'paired_p_r_e'}, {'paired_p_o_s_t'}, {'unpaired_p_r_e'}, {'unpaired_p_o_s_t'}], 1, [0.35, 0.35, 0.35]);
    
    ylabel('single pulse resp (mean dF/F)')
    fig_wrapup(8, script_name);
    
    
    %computing pre-post difference in peak response for paired and unpaired odors
    pk_diff_simple = cat(2, (saved_resps_simple(:, 1) - saved_resps_simple(:, 2)), (saved_resps_simple(:, 3) - saved_resps_simple(:, 4)));
    pk_diff_transition = cat(2, (saved_resps_hover(:, 1) - saved_resps_hover(:, 2)), (saved_resps_hover(:, 3) - saved_resps_hover(:, 4)));
    
    marker_colors = [[0.6, 0.6, 0.6]; [0.6, 0.6, 0.6]];
    col_pairs = [1, 2];
    scattered_dot_plot(pk_diff_simple, 11, 1, 4, 8, marker_colors, 1, col_pairs, [0.75, 0.75, 0.75],...
                            [{'paired'}, {'unpaired'}], 1, [0.9, 0.3, 0.3]);
    ylabel('pre-post peak difference, simple (mean dF/F)')
    fig_wrapup(11, script_name);
    
    
    marker_colors = [[0.6, 0.6, 0.6]; [0.6, 0.6, 0.6]];
    col_pairs = [1, 2];
    scattered_dot_plot(pk_diff_transition, 12, 1, 4, 8, marker_colors, 1, col_pairs, [0.75, 0.75, 0.75],...
                            [{'paired'}, {'unpaired'}], 1, [0.9, 0.3, 0.3]);
    ylabel('pre-post peak difference, transition (mean dF/F)')
    fig_wrapup(12, script_name);
    
    
    %statistical testing
    [h_paired_s, p_paired_s] = ttest(saved_resps_simple(:, 1), saved_resps_simple(:, 2))
    [h_unpaired_s, p_unpaired_s] = ttest(saved_resps_simple(:, 3), saved_resps_simple(:, 4))
    
    [h_paired_t, p_paired_t] = ttest(saved_resps_hover(:, 1), saved_resps_hover(:, 2))
    [h_unpaired_t, p_unpaired_t] = ttest(saved_resps_hover(:, 3), saved_resps_hover(:, 4))
    
    [h_pkdiff_s, p_pkdiff_s] = ttest(pk_diff_simple(:, 1), pk_diff_simple(:, 2))
    [h_pkdiff_t, p_pkdiff_t] = ttest(pk_diff_transition(:, 1), pk_diff_transition(:, 2))
    
    
    %comparing paired odor with unpaired odor
    [h_pre_comp, p_pre_comp] = ttest(saved_resps_hover(:, 1), saved_resps_hover(:, 3))
    [h_post_comp, p_post_comp] = ttest(saved_resps_hover(:, 2), saved_resps_hover(:, 4))
    
    
    
    %computing area of difference between pre and post response traces and sorting flies
    diff_mat_simple  = squeeze(mean(cat(2, (saved_mean_traces_simple(:,1, :) - saved_mean_traces_simple(:,2, :)),  (saved_mean_traces_simple(:,3, :) - saved_mean_traces_simple(:,4, :))), 1));
    diff_mat_transition  = squeeze(mean(cat(2, (saved_mean_traces_transition(:,1, :) - saved_mean_traces_transition(:,2, :)),  (saved_mean_traces_transition(:,3, :) - saved_mean_traces_transition(:,4, :))), 1));
   
    marker_colors = [[0.6, 0.6, 0.6]; [0.6, 0.6, 0.6]];
    col_pairs = [1, 2];
    scattered_dot_plot(diff_mat_simple', 9, 1, 4, 8, marker_colors, 1, col_pairs, [0.75, 0.75, 0.75],...
                            [{'paired'}, {'unpaired'}], 1, [0.9, 0.3, 0.3]);
    ylabel('pre-post area difference, simple (mean dF/F)')
    fig_wrapup(9, script_name);
    
    
    marker_colors = [[0.6, 0.6, 0.6]; [0.6, 0.6, 0.6]];
    col_pairs = [1, 2];
    scattered_dot_plot(diff_mat_transition', 10, 1, 4, 8, marker_colors, 1, col_pairs, [0.75, 0.75, 0.75],...
                            [{'paired'}, {'unpaired'}], 1, [0.9, 0.3, 0.3]);
    ylabel('pre-post area difference, transition (mean dF/F)')
    fig_wrapup(10, script_name);
    
    

end