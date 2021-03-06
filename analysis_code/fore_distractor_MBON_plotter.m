clear all
close all

dataset_list_paths = [...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONAlpha1.xls'};...
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

fore_colour = [0.5, 0.4, 1];
distr_colour = [0.5, 1, 0.4];

for list_n = 1:size(dataset_list_paths, 1)
    curr_dir_list_path = dataset_list_paths{list_n, 1};
    [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
    n_dirs = size(dir_list, 1);
    dataset_list_name = findstr(curr_dir_list_path, 'list_');
    dataset_list_name = curr_dir_list_path((dataset_list_name + 5):(end - 4));
    
    saved_resps = zeros(n_dirs, 5) + nan;                    %mean response of each fly
    
    %loop to go through all experiment datasets listed in list file
    for dir_n = 1:n_dirs
        fly_n = fly_n + 1;
              
        saved_an_results.scriptname = mfilename('fullpath');
        curr_dir = [dir_list{dir_n, 1}, '\'];
        curr_dir = manage_base_paths(curr_dir, 2);
       
        tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
        tif_times = tif_times.time_stamps;
        [stim_mat, stim_mat_simple, column_heads, color_vec] = load_params_trains_modular(curr_dir, tif_times);    %reading in trial stimulus parameters after matching time stamps to F traces
        odor_colors = [color_vec(3, :); color_vec(3, :).*0.75; color_vec(2, :)];
       
        %Reading in experimental parameters
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
        del = find(dff_data_mat_f < -1);
        dff_data_mat_f(del) = -1;       %forcing crazy values to sane ones
        
        %identifying lists of trial numbers
        led_on_col_n = find_stim_mat_simple_col('led_on', column_heads);            %identifying relevant column number in stim_mat_simple
        od_olf1_col_n = find_stim_mat_simple_col('odor_n', column_heads);           %identifying relevant column number in stim_mat_simple
        od_olf2_col_n = find_stim_mat_simple_col('odour_olf2', column_heads);       %identifying relevant column number in stim_mat_simple
        dur_olf1_col_n = find_stim_mat_simple_col('duration', column_heads);        %identifying relevant column number in stim_mat_simple
        dur_olf2_col_n = find_stim_mat_simple_col('duration_olf2', column_heads);   %identifying relevant column number in stim_mat_simple
        od_col_ns = [od_olf1_col_n, od_olf2_col_n];
        dur_col_ns = [dur_olf1_col_n, dur_olf2_col_n];
        
        pairing_tr_n = find(stim_mat_simple(:, led_on_col_n) == 1);
        [stim_frs_pairing_tr, fore_olf_n] = compute_stim_frs_modular(stim_mat, pairing_tr_n, frame_time);        %olfactometer number used to deliver foreground odor in pairing trial
        distr_olf_n = [1, 2];
        distr_olf_n(distr_olf_n == fore_olf_n) = [];                 %olfactometer number used to deliver distractor odor in pairing trial
        saved_resps(dir_n, 5) = fore_olf_n;
        
        %assigning column numbers in stim_mat_simple depedning on the
        %olfactometer used to delvier the foreground odour in the pairing
        %trial.
        if fore_olf_n == 1
            fore_od_col_n = od_olf1_col_n;
            fore_dur_col_n = dur_olf1_col_n;
            distr_od_col_n = od_olf2_col_n;
            distr_dur_col_n = dur_olf2_col_n;
        elseif fore_olf_n == 2
            fore_od_col_n = od_olf2_col_n;
            fore_dur_col_n = dur_olf2_col_n;
            distr_od_col_n = od_olf1_col_n;
            distr_dur_col_n = dur_olf1_col_n;            
        else
        end
        
        fore_od_n = stim_mat_simple(pairing_tr_n, od_col_ns(fore_olf_n) );
        distr_od_n = stim_mat_simple(pairing_tr_n, od_col_ns(distr_olf_n) );
                
        fore_od_trs_temp = find(stim_mat_simple(:, fore_od_col_n) == fore_od_n);        %preliminary list of foreground odor presentation trials
        fore_od_trs_temp = fore_od_trs_temp(fore_od_trs_temp < pairing_tr_n);
        pre_post_dur = stim_mat_simple(fore_od_trs_temp, fore_dur_col_n);               %list of durations that foreground odor was delivered on
        pre_post_dur = max(pre_post_dur, [], 'omitnan');                                    %getting rid of very small durations used as dummies for olfactometer1.
        
        fore_od_trs = find(stim_mat_simple(:, fore_od_col_n) == fore_od_n & stim_mat_simple(:, fore_dur_col_n) == pre_post_dur);             %list of foreground odor presentation trials
        distr_od_trs = find(stim_mat_simple(:, distr_od_col_n) == distr_od_n & stim_mat_simple(:, distr_dur_col_n) == pre_post_dur);         %list of foreground odor presentation trials
        
        
        %1. pre, foreground odor resps
        curr_trs = fore_od_trs(fore_od_trs < pairing_tr_n);
        curr_trs = sort(curr_trs);
        curr_trs(1) = [];             %getting rid of first preesntation of foreground odour     
        curr_traces = dff_data_mat(:, :, curr_trs);
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
        stim_frs = stim_frs{fore_olf_n};
        saved_resps(dir_n, 1) = mean(mean(squeeze(curr_traces(stim_frs(1):(stim_frs(2) + round(4./frame_time)), :, :)), 1, 'omitnan'));
        trace_mean = mean(curr_traces, 3, 'omitnan');
        trace_se = std(curr_traces, [], 3, 'omitnan')./sqrt(length(curr_trs));
        
        figure(1)
%         plot(trace_mean, 'Color', fore_colour, 'lineWidth', 3)
%         hold on
%         plot(squeeze(curr_traces), 'Color', fore_colour, 'lineWidth', 0.2)
        shadedErrorBar([], trace_mean, trace_se, {'Color', fore_colour})
        hold on
        ylabel('foreground odor responses (dF/F)')
        set_xlabels_time(1, frame_time, 10);
        
        %2. pre, distractor odor resps
        curr_trs = distr_od_trs(distr_od_trs < pairing_tr_n);
        curr_trs = sort(curr_trs);
        curr_trs(1) = [];             %getting rid of first preesntation of foreground odour     
        curr_traces = dff_data_mat(:, :, curr_trs);
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
        stim_frs = stim_frs{fore_olf_n};
        saved_resps(dir_n, 2) = mean(mean(squeeze(curr_traces(stim_frs(1):(stim_frs(2) + round(4./frame_time)), :, :)), 1, 'omitnan'));
        trace_mean = mean(curr_traces, 3, 'omitnan');
        trace_se = std(curr_traces, [], 3, 'omitnan')./sqrt(length(curr_trs));
        
        figure(2)
%         plot(trace_mean, 'Color', distr_colour, 'lineWidth', 3)
%         hold on
%         plot(squeeze(curr_traces), 'Color', distr_colour, 'lineWidth', 0.2)
        shadedErrorBar([], trace_mean, trace_se, {'Color', distr_colour})
        hold on
        ylabel('distractor odor responses (dF/F)')
        set_xlabels_time(2, frame_time, 10);
        
        
        %3. post, foreground odor resps
        curr_trs = fore_od_trs(fore_od_trs > pairing_tr_n);
        curr_trs = sort(curr_trs);
        curr_traces = dff_data_mat(:, :, curr_trs);
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
        stim_frs = stim_frs{fore_olf_n};
        saved_resps(dir_n, 3) = mean(mean(squeeze(curr_traces(stim_frs(1):(stim_frs(2) + round(4./frame_time)), :, :)), 1, 'omitnan'));
        trace_mean = mean(curr_traces, 3, 'omitnan');
        trace_se = std(curr_traces, [], 3, 'omitnan')./sqrt(length(curr_trs));
        
        figure(1)
%         plot(trace_mean, 'Color', fore_colour.*0.7, 'lineWidth', 3)
%         hold on
%         plot(squeeze(curr_traces), 'Color', fore_colour.*0.7, 'lineWidth', 0.2)
        shadedErrorBar([], trace_mean, trace_se, {'Color', fore_colour.*0.7})
        fig_wrapup(1, script_name);
        add_stim_bar(1, stim_frs, [0.75, 0.75, 0.75]);
        hold off        
        
        %4. post, distractor odor resps
        curr_trs = distr_od_trs(distr_od_trs > pairing_tr_n);
        curr_trs = sort(curr_trs);
        curr_traces = dff_data_mat(:, :, curr_trs);
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
        stim_frs = stim_frs{fore_olf_n};
        saved_resps(dir_n, 4) = mean(mean(squeeze(curr_traces(stim_frs(1):(stim_frs(2) + round(4./frame_time)), :, :)), 1, 'omitnan'));
        trace_mean = mean(curr_traces, 3, 'omitnan');
        trace_se = std(curr_traces, [], 3, 'omitnan')./sqrt(length(curr_trs));
        
        figure(2)
%         plot(trace_mean, 'Color', distr_colour.*0.7, 'lineWidth', 3)
%         hold on
%         plot(squeeze(curr_traces), 'Color', distr_colour.*0.7, 'lineWidth', 0.2)
        shadedErrorBar([], trace_mean, trace_se, {'Color', distr_colour.*0.7})
        fig_wrapup(2, script_name);
        add_stim_bar(2, stim_frs, [0.75, 0.75, 0.75]);
        hold off
        
        if suppress_plots == 0
            keyboard
        else
        end
        
    end
    
    keyboard
end