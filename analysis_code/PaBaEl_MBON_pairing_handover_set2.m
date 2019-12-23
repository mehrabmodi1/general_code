clear all
close all

dataset_list_paths = [...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONAlpha1.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONAlpha1_set2.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONAlpha1_set3_highLED.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONGamma2.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONGamma2_set1.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONGamma2_set2.xls'};...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_PaBaEl_MBONG2_handover.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_MBONG2_PaBaEl_handover_simple_starved.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_MBONG2_PaBaEl_handover_starved_set2.xls'};...
                        
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

for list_n = 1:size(dataset_list_paths, 1)
    curr_dir_list_path = dataset_list_paths{list_n, 1};
    [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
    n_dirs = size(dir_list, 1);
    dataset_list_name = findstr(curr_dir_list_path, 'list_');
    dataset_list_name = curr_dir_list_path((dataset_list_name + 5):(end - 4));
    
    saved_resps = zeros(n_dirs, 7) + nan;                    %mean response of each fly
    
    %loop to go through all experiment datasets listed in list file
    saved_resps = [];
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
        raw_data_mat = raw_data_mat(:, :, stim_mat_simple(:, tif_n_col_n));       %making sure only time-stamp matched trials are used for further analysis
        n_cells = size(raw_data_mat, 2);
        
        %calculating dF/F traces from raw data
        filt_time = 0.5;            %in s, the time window for boxcar filter for generating filtered traces
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
        
        unpaired_od_n_olf2 = odn_list_olf2;
        unpaired_od_n_olf2(unpaired_od_n_olf2 == intersect(odn_list_olf2, paired_od_n_olf2)) = [];
        unpaired_od_name = odor_names2{unpaired_od_n_olf2};
        unpaired_od_n_olf1 = od_name_lookup(odor_names1, unpaired_od_name); 
        
        handover_trs = find(stim_mat_simple(:, dur_col_ns(1)) == od_durs(1) & stim_mat_simple(:, dur_col_ns(2)) == od_durs(1));
        
        %building tr lists
        all_trs_paired_od_olf2 = find(stim_mat_simple(:, od_col_ns(2)) == paired_od_n_olf2 & stim_mat_simple(:, dur_col_ns(2)) == od_durs(1));   %this includes post-pairing and handover trials
        all_trs_unpaired_od_olf1 = find(stim_mat_simple(:, od_col_ns(1)) == unpaired_od_n_olf1 & stim_mat_simple(:, dur_col_ns(1)) == od_durs(1));   %this includes post-pairing and handover trials
        
        %1. non-handover, paired odor
        curr_trs = all_trs_paired_od_olf2;
        [del, del2] = intersect(curr_trs, handover_trs);
        curr_trs(del2) = [];            %excluded handover trs
        curr_trs_pre = curr_trs(curr_trs < pairing_tr_n);   %only pre trs 
        curr_trs_post = curr_trs(curr_trs > pairing_tr_n);   %only post trs 
                
        curr_traces_pre = dff_data_mat_f(:, :, curr_trs_pre);
        curr_traces_post = dff_data_mat_f(:, :, curr_trs_post);
        paired_pre_mean = median(squeeze(curr_traces_pre), 2, 'omitnan');
        paired_post_mean = median(squeeze(curr_traces_post), 2, 'omitnan');
        length_data_paired = length(paired_post_mean) - sum(isnan(paired_post_mean));
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs_pre(1), frame_time);
        stim_frs = stim_frs{2};                   %because olf2 is used for paired odor
        
        figure(1)
        plot(squeeze(curr_traces_pre), 'Color', [0.6, 0.6, 0.6], 'lineWidth', 3);
        hold on
        plot(squeeze(curr_traces_post), 'Color', [0.1, 0.1, 0.1], 'lineWidth', 3);
        %shadedErrorBar([], trace_mean, trace_se, {'Color', fore_colour})
        ylabel('paired odor, simple (dF/F)')
        ax_vals = axis;
        ax_vals(1, 4) = 0.1;
        axis(ax_vals);
        set_xlabels_time(1, frame_time, 10);
        fig_wrapup(1, script_name);
        add_stim_bar(1, stim_frs, color_vec(1, :));
        hold off
        
        
        %3. non-handover, unpaired odor
        curr_trs = all_trs_unpaired_od_olf1;
        [del, del2] = intersect(curr_trs, handover_trs);
        curr_trs(del2) = [];            %excluded handover trs
        curr_trs_pre = curr_trs(curr_trs < pairing_tr_n);   %only pre trs 
        curr_trs_post = curr_trs(curr_trs > pairing_tr_n);   %only post trs 
                
        curr_traces_pre = dff_data_mat_f(:, :, curr_trs_pre);
        curr_traces_post = dff_data_mat_f(:, :, curr_trs_post);
        unpaired_pre_mean = median(squeeze(curr_traces_pre), 2);
        unpaired_post_mean = median(squeeze(curr_traces_pre), 2);
        length_data_unpaired = length(unpaired_post_mean) - sum(isnan(unpaired_post_mean));
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs_pre(1), frame_time);
        stim_frs = stim_frs{1};                   %because olf2 is used for paired odor
        
        
        figure(3)
        plot(squeeze(curr_traces_pre), 'Color', [0.6, 0.6, 0.6], 'lineWidth', 3);
        hold on
        plot(squeeze(curr_traces_post), 'Color', [0.1, 0.1, 0.1], 'lineWidth', 3);
        ylabel('unpaired odor, simple (dF/F)')
        ax_vals = axis;
        ax_vals(1, 4) = 0.1;
        axis(ax_vals);
        set_xlabels_time(3, frame_time, 10);
        fig_wrapup(3, script_name);
        add_stim_bar(3, stim_frs, color_vec(2, :));
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
        curr_trs_pre = curr_trs(curr_trs < pairing_tr_n);   %only pre trs 
        curr_trs_post = curr_trs(curr_trs > pairing_tr_n);   %only post trs 
                
        curr_traces_pre = dff_data_mat_f(:, :, curr_trs_pre);
        curr_traces_post = dff_data_mat_f(:, :, curr_trs_post);
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs_pre(1), frame_time);
        stim_frs = [stim_frs{1}; stim_frs{2}];                   %because olf2 is used for paired odor
        
        %computing response to second odor pulse in handover presentations
        paired_mean_resp_pre = mean(mean(squeeze(curr_traces_pre(stim_frs(2, 1):stim_frs(2, 2), :)), 2));
        paired_mean_resp_post = mean(mean(squeeze(curr_traces_post(stim_frs(2, 1):stim_frs(2, 2), :)), 2));
        
        %raw trace plot
        figure(2)
        plot(squeeze(curr_traces_pre), 'Color', [0.6, 0.6, 0.6], 'lineWidth', 3);
        hold on
        plot(squeeze(curr_traces_post), 'Color', [0.1, 0.1, 0.1], 'lineWidth', 3);
        ylabel('paired odor, handover (dF/F)')
        ax_vals = axis;
        ax_vals(1, 4) = 0.1;
        axis(ax_vals);
        set_xlabels_time(2, frame_time, 10);
        fig_wrapup(2, script_name);
        add_stim_bar(2, stim_frs, [color_vec(2, :); color_vec(1, :);]);
        hold off
        
        
        %linear model plot
        figure(5)
        mean_pre = mean(squeeze(curr_traces_pre), 2, 'omitnan');
        se_pre = std(squeeze(curr_traces_pre), [], 2, 'omitnan')./sqrt(size(squeeze(curr_traces_pre), 2));
        mean_post = mean(squeeze(curr_traces_post), 2, 'omitnan');
        se_post = std(squeeze(curr_traces_post), [], 2, 'omitnan')./sqrt(size(squeeze(curr_traces_post), 2));
        shadedErrorBar([], mean_pre, se_pre, {'Color', [0.6, 0.6, 0.6], 'lineWidth', 1.5});
        hold on
        shadedErrorBar([], mean_post, se_post, {'Color', [0, 0, 0], 'lineWidth', 1.5});
        plot(minusplus_sum_pre, 'Color', [0.5, 0.7, 0.9], 'lineWidth', 1.5);
        plot(minusplus_sum_post, 'Color', [0.2, 0.35, 0.6], 'lineWidth', 1.5);
        ylabel('paired odor, handover (dF/F)')
        ax_vals = axis;
        ax_vals(1, 4) = 0.1;
        axis(ax_vals);
        set_xlabels_time(5, frame_time, 10);
        fig_wrapup(5, script_name);
        add_stim_bar(5, stim_frs, [color_vec(2, :); color_vec(1, :);]);
        hold off
        
        %4. handover, unpaired odor
        curr_trs = handover_trs;
        [del, del2] = intersect(curr_trs, all_trs_unpaired_od_olf1);
        curr_trs(del2) = [];            %got rid of handover trs with unpaired od delivered first on olf1
        curr_trs_pre = curr_trs(curr_trs < pairing_tr_n);   %only pre trs 
        curr_trs_post = curr_trs(curr_trs > pairing_tr_n);   %only post trs 
                
        curr_traces_pre = dff_data_mat_f(:, :, curr_trs_pre);
        curr_traces_post = dff_data_mat_f(:, :, curr_trs_post);
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs_pre(1), frame_time);
        stim_frs = [stim_frs{1}; stim_frs{2}];                   %because olf2 is used for paired odor
        
        %computing response to second odor pulse in handover presentations
        unpaired_mean_resp_pre = mean(mean(squeeze(curr_traces_pre(stim_frs(2, 1):stim_frs(2, 2), :)), 2));
        unpaired_mean_resp_post = mean(mean(squeeze(curr_traces_post(stim_frs(2, 1):stim_frs(2, 2), :)), 2));
        
        saved_resps = [saved_resps; paired_mean_resp_pre, paired_mean_resp_post, unpaired_mean_resp_pre, unpaired_mean_resp_post];
        
        figure(4)
        plot(squeeze(curr_traces_pre), 'Color', [0.6, 0.6, 0.6], 'lineWidth', 3);
        hold on
        plot(squeeze(curr_traces_post), 'Color', [0.1, 0.1, 0.1], 'lineWidth', 3);
        ylabel('unpaired odor, handover (dF/F)')
        ax_vals = axis;
        ax_vals(1, 4) = 0.1;
        axis(ax_vals);
        set_xlabels_time(4, frame_time, 10);
        fig_wrapup(4, script_name);
        add_stim_bar(4, stim_frs, [color_vec(1, :); color_vec(2, :);]);
        hold off
        
        
        %linear model plot
        figure(6)
        mean_pre = mean(squeeze(curr_traces_pre), 2, 'omitnan');
        se_pre = std(squeeze(curr_traces_pre), [], 2, 'omitnan')./sqrt(size(squeeze(curr_traces_pre), 2));
        mean_post = mean(squeeze(curr_traces_post), 2, 'omitnan');
        se_post = std(squeeze(curr_traces_post), [], 2, 'omitnan')./sqrt(size(squeeze(curr_traces_post), 2));
        shadedErrorBar([], mean_pre, se_pre, {'Color', [0.6, 0.6, 0.6], 'lineWidth', 1.5});
        hold on
        shadedErrorBar([], mean_post, se_post, {'Color', [0, 0, 0], 'lineWidth', 1.5});
        plot(plusminus_sum_pre, 'Color', [0.5, 0.7, 0.9], 'lineWidth', 1.5);
        plot(plusminus_sum_post, 'Color', [0.2, 0.35, 0.6], 'lineWidth', 1.5);
        ylabel('unpaired odor, handover (dF/F)')
        ax_vals = axis;
        ax_vals(1, 4) = 0.1;
        axis(ax_vals);
        set_xlabels_time(6, frame_time, 10);
        fig_wrapup(6, script_name);
        add_stim_bar(6, stim_frs, [color_vec(1, :); color_vec(2, :);]);
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
    marker_colors = [color_vec(1, :); color_vec(1, :); color_vec(2, :); color_vec(2, :)];
    col_pairs = [1, 2; 3, 4];
    scattered_dot_plot(saved_resps(:, 1:4), 7, 1, 4, 8, marker_colors, 1, col_pairs, [0.75, 0.75, 0.75],...
                            [{'paired_p_r_e'}, {'paired_p_o_s_t'}, {'unpaired_p_r_e'}, {'unpaired_p_o_s_t'}], 1, [0.35, 0.35, 0.35]);
    
    ylabel('response size (mean dF/F)')
    fig_wrapup(7, script_name);
    
    %statistical testing
    [h_paired, p_paired] = ttest(saved_resps(:, 1), saved_resps(:, 2));
    [h_unpaired, p_unpaired] = ttest(saved_resps(:, 3), saved_resps(:, 4));
    
    keyboard
    
    

end