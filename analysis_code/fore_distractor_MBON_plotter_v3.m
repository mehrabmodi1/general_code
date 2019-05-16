clear all
close all

dataset_list_paths = [...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONAlpha1.xls'};...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_fore_distr_MBONAlpha1_set2.xls'};...
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
    
    saved_resps = zeros(n_dirs, 7) + nan;                    %mean response of each fly
    
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
        
        %identifying stim_mat_simple col numbers
        led_on_col_n = find_stim_mat_simple_col('led_on', column_heads);            %identifying relevant column number in stim_mat_simple
        od_olf1_col_n = find_stim_mat_simple_col('odor_n', column_heads);           %identifying relevant column number in stim_mat_simple
        od_olf2_col_n = find_stim_mat_simple_col('odour_olf2', column_heads);       %identifying relevant column number in stim_mat_simple
        dur_olf1_col_n = find_stim_mat_simple_col('duration', column_heads);        %identifying relevant column number in stim_mat_simple
        dur_olf2_col_n = find_stim_mat_simple_col('duration_olf2', column_heads);   %identifying relevant column number in stim_mat_simple
        od_col_ns = [od_olf1_col_n, od_olf2_col_n];
        dur_col_ns = [dur_olf1_col_n, dur_olf2_col_n];
        
        pairing_tr_n = find(stim_mat_simple(:, led_on_col_n) == 1);
        fore_od_n_olf1 = stim_mat_simple(pairing_tr_n, od_olf1_col_n);          %This assumes the foreground odor is always delivered on olf1
        distr_od_n_olf2 = stim_mat_simple(pairing_tr_n, od_olf2_col_n);
        saved_resps(dir_n, 7) = fore_od_n_olf1;                                 %keeping track of the foreground odor number for current dataset
        
        %identifying matching odor_ns for olf1 and olf2 based on odorNames
        %for the two olfactometers
        fore_od_name = odor_names1{fore_od_n_olf1};
        fore_od_n_olf2 = lookup_cell_vec(fore_od_name, odor_names2);
        distr_od_name = odor_names2{distr_od_n_olf2};
        distr_od_n_olf1 = lookup_cell_vec(distr_od_name, odor_names1);
        fore_od_n = fore_od_n_olf1;         %olf1 is always used to deliver the foreground odor
        distr_od_n = distr_od_n_olf2;        
        od_list_olf1 = unique(stim_mat_simple(:, od_olf1_col_n));                    %all odors delivered on olf1
        ctrl_od_ni = find(od_list_olf1 ~= fore_od_n_olf1 & od_list_olf1 ~= distr_od_n_olf1);
        ctrl_od_n_olf1 = od_list_olf1(ctrl_od_ni);
        
        %checking which olfactometer was used to deliver the foreground
        %odor in the pre-post trials.
        fore_trs_olf2 = find(stim_mat_simple(1:(pairing_tr_n - 1), od_olf2_col_n) == fore_od_n_olf2);
        if isempty(fore_trs_olf2) ==  1
            fore_olf_n = 1;                 %olfactometer1 used for foreground odor delivery in pre-post trials
            distr_olf_n = 2;                %olfactometer2 used for distractor odor delivery in pre-post trials
            fore_od_col_n = od_olf1_col_n;
            fore_dur_col_n = dur_olf1_col_n;
            distr_od_col_n = od_olf2_col_n;
            distr_dur_col_n = dur_olf2_col_n;            
        else
            fore_olf_n = 2;                 %olfactometer2 used for foreground odor delivery in pre-post trials
            distr_olf_n = 1;                %olfactometer1 used for distractor odor delivery in pre-post trials
            fore_od_col_n = od_olf2_col_n;
            fore_dur_col_n = dur_olf2_col_n;
            distr_od_col_n = od_olf1_col_n;
            distr_dur_col_n = dur_olf1_col_n;
            
        end
        
        pre_post_dur = stim_mat_simple(1:(pairing_tr_n - 1), fore_dur_col_n);               
        pre_post_dur = max(pre_post_dur, [], 'omitnan');                        %duration that pre-post odors were delivered for
        
        fore_od_trs = find(stim_mat_simple(:, fore_od_col_n) == fore_od_n & stim_mat_simple(:, fore_dur_col_n) == pre_post_dur);             %list of foreground odor presentation trials
        distr_od_trs = find(stim_mat_simple(:, distr_od_col_n) == distr_od_n & stim_mat_simple(:, distr_dur_col_n) == pre_post_dur);         %list of foreground odor presentation trials
        ctrl_od_trs = find(stim_mat_simple(:, od_olf1_col_n) == ctrl_od_n_olf1);
        
        
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
        stim_frs = stim_frs{distr_olf_n};
        saved_resps(dir_n, 3) = mean(mean(squeeze(curr_traces(stim_frs(1):(stim_frs(2) + round(4./frame_time)), :, :)), 1, 'omitnan'));
        
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
        saved_resps(dir_n, 2) = mean(mean(squeeze(curr_traces(stim_frs(1):(stim_frs(2) + round(4./frame_time)), :, :)), 1, 'omitnan'));
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
        stim_frs = stim_frs{distr_olf_n};
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
        close figure 1
        close figure 2
        
        %5 Ctrl odor pre responses
        curr_trs = ctrl_od_trs(ctrl_od_trs < pairing_tr_n);
        curr_trs = sort(curr_trs);
        curr_traces = dff_data_mat(:, :, curr_trs);
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
        stim_frs = stim_frs{fore_olf_n};
        saved_resps(dir_n, 5) = mean(mean(squeeze(curr_traces(stim_frs(1):(stim_frs(2) + round(4./frame_time)), :, :)), 1, 'omitnan'));
        
        %6 Ctrl odor post responses
        curr_trs = ctrl_od_trs_od_trs(ctrl_od_trs > pairing_tr_n);
        curr_trs = sort(curr_trs);
        curr_traces = dff_data_mat(:, :, curr_trs);
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
        stim_frs = stim_frs{fore_olf_n};
        saved_resps(dir_n, 6) = mean(mean(squeeze(curr_traces(stim_frs(1):(stim_frs(2) + round(4./frame_time)), :, :)), 1, 'omitnan'));
        
        
    end
    marker_colors = [fore_colour; fore_colour; distr_colour; distr_colour; ctrl_colour; ctrl_colour];
    col_pairs = [1, 2; 3, 4; 5, 6];
    scattered_dot_plot(saved_resps(:, 1:6), 5, 1, 4, 8, marker_colors, 1, col_pairs, [0.75, 0.75, 0.75],...
                            [{'fore pre'}, {'fore post'}, {'distr pre'}, {'distr post'}, {'ctrl_pre'}, {'ctrl_post'}], 1, [0.35, 0.35, 0.35]);
    
    fig_wrapup(5, script_name);
    keyboard
end