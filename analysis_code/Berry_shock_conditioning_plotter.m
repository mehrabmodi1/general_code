clear all
close all

dataset_list_paths = [...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\Berry_handover_WTalpha3_fedshock_BAplus.xls'};...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\Berry_handover_WTalpha3_fedshock_PAplus.xls'};...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\Berry_handover_WTalpha3_fedshock_unpairedctrl.xls'};...
                      ];
            
suppress_plots = 1;
plotting_quant_no_filt = 0;     %1 - only unfiltered traces used for all analysis and plotting - traces included. 0 - filtered traces used for everything.
cell_n = 4;

% [del, odor_names1] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
% [del, odor_names2] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList_olf2.xls', 1);
% odor_names2{3} = 'Butyl acetate';

global color_vec;                
a = colormap('bone');
global greymap
greymap = flipud(a);
fly_n = 0;
script_name = mfilename;

y_ax_traces = 0.8;
y_ax_fit_traces = 0.6;

an_save_path = 'C:\Data\Data\Analysed_data\Analysis_results\PaBaEl_Gamma2\';
force_resave = 1;

n_sec = 2;      %width of moving integration window in s

integ_win_s = 5;        %width of pulse integration window for response quantification

for list_n = 1:size(dataset_list_paths, 1)
    
    curr_dir_list_path = dataset_list_paths{list_n, 1};
    [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
    n_dirs = size(dir_list, 1);
    dataset_list_name = findstr(curr_dir_list_path, 'list_');
    dataset_list_name = curr_dir_list_path((dataset_list_name + 5):(end - 4));
    
    if isempty(findstr(curr_dir_list_path, 'BA')) == 0
        paired_od_n_olf2 = 3;
        unpaired_od_n_olf2 = 1;
        set_list_type = 1;      %paired list
    elseif isempty(findstr(curr_dir_list_path, 'PA')) == 0
        paired_od_n_olf2 = 1;
        unpaired_od_n_olf2 = 3;
        set_list_type = 1;      %paired list
    elseif  isempty(findstr(curr_dir_list_path, 'unpaired')) == 0
        set_list_type = 0;  %unpaired list
    end
    
    
    saved_traces_all = [];
    saved_PID_traces_all = [];
    %loop to go through all experiment datasets listed in list file
    for dir_n = 1:n_dirs
        fly_n = fly_n + 1;
             
        curr_dir = [dir_list{dir_n, 1}, '\'];
        curr_dir = manage_base_paths(curr_dir, 2)
       
        tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
        tif_times = tif_times.time_stamps;
        
        cd(curr_dir);
        tif_name = dir('*.tif');
        if isempty(tif_name) == 1
            continue
        else
        end
        stack_obj = ScanImageTiffReader([curr_dir, tif_name(1).name]);
        [frame_time, zoom, n_chans, PMT_offsets] = SI_tif_info(stack_obj);
        
        [stim_mat, stim_mat_simple, column_heads, color_vec, good_tr_list, params_orig, PID_traces] = load_params_trains_modular(curr_dir, tif_times, frame_time);    %reading in trial stimulus parameters after matching time stamps to F traces
        
        %looking up odor numbers across olfactometers based on current dataset's odor name list
        odor_names1 = stim_mat.odourNames;
        odor_names2 = stim_mat.odourNames_olf2;
        paired_od_n_olf1 = od_name_lookup(odor_names1, odor_names2{paired_od_n_olf2});
        unpaired_od_n_olf1 = od_name_lookup(odor_names1, odor_names2{unpaired_od_n_olf2});
        
        paired_color = color_vec(2, :);
        unpaired_color = color_vec(1, :);
        EL_color = color_vec(3, :);
        %mean_color = [0.5, 0.83, 0.98];
        %mean_color = [149, 200, 216]./256;
        mean_color = ([0, 49, 152]./256).*1.5;
        
        stim_mat_simple_nonans = stim_mat_simple;
        stim_mat_simple_nonans(isnan(stim_mat_simple)) = 0;
                
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
                
        %Reading in experimental parameters
        odor_list_olf1 = unique(stim_mat_simple(:, od_col_ns(1) ) );
        n_odors_olf1 = length(odor_list_olf1);
        odor_dur_list_olf1 = unique(stim_mat_simple(:, dur_col_ns(1) ) );
        n_od_durs_olf1 = length(odor_dur_list_olf1);
                
        odor_list_olf2 = unique(stim_mat_simple(:, od_col_ns(2) ) );
        n_odors_olf2 = length(odor_list_olf2);
        odor_dur_list_olf2 = unique(stim_mat_simple(:, dur_col_ns(2) ) );
        n_od_durs_olf2 = length(odor_dur_list_olf2);
        
        
        integ_win = round(integ_win_s./frame_time); %integration window for response quantification in frames
       
        %loading extracted raw fluorescence data matrices written by raw_dff_extractor
        raw_data_mat = load([curr_dir 'extracted_raw_data_mat.mat']);
        raw_data_mat = raw_data_mat.raw_data_mat;           %raw F traces extracted from ROIs
        raw_data_mat_orig = raw_data_mat;
        tif_n_col_n = find_stim_mat_simple_col('matched_tif_n', column_heads);
        n_cells = size(raw_data_mat, 2);
        if cell_n > n_cells
            cell_n = 1;
            disp('WARNING! specified cell_n not valid, reassigning cell_n = 1.');
        else
        end
        
        %inserting dummy trials (all nans) into raw_data_mat for pairing trials for
        %which no corress .tifs were acquired
        [raw_data_mat, good_tr_list] = match_up_rawmat_matchedtrs(raw_data_mat, stim_mat_simple, tif_n_col_n, good_tr_list);
        
        %if weird EL trial1 exists, getting rid of it
        if stim_mat_simple(1, od_olf1_col_n) == 11
            stim_mat(1) = [];
            stim_mat_simple(1, :) = [];
            raw_data_mat(:, :, 1) = [];
        else
        end
        
        n_cells = size(raw_data_mat, 2);
        
        %calculating dF/F traces from raw data
        filt_time = 0.2;            %in s, the time window for boxcar filter for generating filtered traces
        [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, curr_dir, good_tr_list);
        if plotting_quant_no_filt == 1
            dff_data_mat_f = dff_data_mat;
        else
        end
          
%         if size(dff_data_mat, 2) > 1
%             dff_data_mat = mean(dff_data_mat, 2, 'omitnan');
%             dff_data_mat_f = mean(dff_data_mat_f, 2, 'omitnan');
%         else
%         end
        dff_data_mat = dff_data_mat(:, cell_n, :);
        dff_data_mat_f  = dff_data_mat_f(:, cell_n, :);
        del = find(dff_data_mat_f < -1);
        dff_data_mat_f(del) = -1;       %forcing crazy values to sane ones
                
        keyboard
                
        y_ax_lim = [];
        plot_means = 1;
        
        %plotting, quantifying resps
        tr_lists = [1:3; 6:8; 11:13];   %cases where EL alone trials exist
               
        %looping through pre, post1 and post2 trial types
        resp_mat = zeros(3, 3);
        
        saved_traces_curr = [];
        for tr_type_n = 1
            tr_list = tr_lists(tr_type_n, :);
            col_multiplier = 0.6.^(tr_type_n - 1);
            
            %paired odor response
            if set_list_type == 0   %unpaired control trials
                curr_tr = find(stim_mat_simple(tr_list, od_olf1_col_n) == paired_od_n_olf1 & stim_mat_simple(tr_list, led_on_col_n) == 0);
                stim_frs_paired = compute_stim_frs_modular(stim_mat, curr_tr(1), frame_time);
                stim_frs = stim_frs_paired{1};      %integrating over entire pulse1 odor period for simple stim protocols
                stim_frs_EL = stim_frs_paired{1};
            elseif set_list_type == 1 
                curr_tr = find(stim_mat_simple(tr_list, od_olf2_col_n) == paired_od_n_olf2 & stim_mat_simple(tr_list, led_on_col_n) == 0);
                stim_frs_paired = compute_stim_frs_modular(stim_mat, curr_tr(1), frame_time);
                stim_frs = stim_frs_paired{2};      %integrating over entire pulse2 odor period for simple stim protocols
                stim_frs_EL = stim_frs_paired{1};
            elseif set_list_type == 3       %for handover, EL second trials
                curr_tr = find(stim_mat_simple(tr_list, od_olf1_col_n) == paired_od_n_olf1 & stim_mat_simple(tr_list, led_on_col_n) == 0);
                stim_frs_paired = compute_stim_frs_modular(stim_mat, curr_tr(1), frame_time);
                stim_frs = stim_frs_paired{2};      %integrating over entire pulse2 odor period for simple stim protocols
                stim_frs_EL = stim_frs_paired{1};
            else
            end
            
            stim_frs(2) = stim_frs(2);  %adding on 2s after odor off to extend integration window
            stim_frs_EL(2) = stim_frs_EL(2); %adding on 2s after odor off to extend integration window
            curr_tr = tr_list(curr_tr);
            
            curr_trace = squeeze(dff_data_mat_f(:, :, curr_tr));
            resp_vec(1, 1) = mean(mean(curr_trace(stim_frs(1):((stim_frs(1) + integ_win) + round(3./frame_time)) )));
            try
                saved_traces_curr(:, 1, tr_type_n) = curr_trace;
            catch
                keyboard
            end
                
                
            saved_PID_traces_curr(:, 1, tr_type_n) = PID_traces(:, curr_tr);
            
            fig_h = figure(1);
            set(fig_h, 'Position', [100, 100, 1200, 350]);
            subplot(1, 3, 1);
            plot(curr_trace, 'lineWidth', 1.5, 'Color', paired_color.*col_multiplier);
            hold on

            %un-paired odor response
            if set_list_type == 0
                curr_tr = find(stim_mat_simple(tr_list, od_olf1_col_n) == unpaired_od_n_olf1);
            elseif set_list_type == 1 || set_list_type == 2     %for handover without EL or handover EL first trials
                curr_tr = find(stim_mat_simple(tr_list, od_olf2_col_n) == unpaired_od_n_olf2);
            elseif set_list_type == 3       %for handover, EL second trials
                curr_tr = find(stim_mat_simple(tr_list, od_olf1_col_n) == unpaired_od_n_olf1);
            else
            end
            
            curr_tr = tr_list(curr_tr);
            curr_trace = squeeze(dff_data_mat_f(:, :, curr_tr));
            resp_vec(1, 2) = mean(mean(curr_trace(stim_frs(1):((stim_frs(1) + integ_win) + round(3./frame_time)) )));
            saved_traces_curr(:, 2, tr_type_n) = curr_trace;
            saved_PID_traces_curr(:, 2, tr_type_n) = PID_traces(:, curr_tr);
            subplot(1, 3, 2)
            plot(curr_trace, 'lineWidth', 1.5, 'Color', unpaired_color.*col_multiplier);
            hold on

            %EL odor response
            if set_list_type <= 1
                curr_tr = find(stim_mat_simple(tr_list, od_olf1_col_n) == 11);
                curr_tr = tr_list(curr_tr);
                curr_trace = squeeze(dff_data_mat_f(:, :, curr_tr));
                resp_vec(1, 3) = mean(mean(curr_trace(stim_frs_EL(1):((stim_frs_EL(1) + integ_win) + round(3./frame_time)) )));
            else
                curr_trace = zeros(size(dff_data_mat_f, 1), 1) + nan;    %since EL alone trials don't exist, padding
                resp_vec(1, 3) = nan; 
            end
            stim_frs_bar_EL = compute_stim_frs_modular(stim_mat, curr_tr(1), frame_time);
            stim_frs_bar_EL = stim_frs_bar_EL{2};
            saved_traces_curr(:, 3, tr_type_n) = curr_trace;
            saved_PID_traces_curr(:, 3, tr_type_n) = PID_traces(:, curr_tr);
            subplot(1, 3, 3)
            plot(curr_trace, 'lineWidth', 1.5, 'Color', EL_color.*col_multiplier);
            hold on
            resp_mat(tr_type_n, :) = resp_vec;
        end
        
        saved_traces_all = pad_n_concatenate(saved_traces_all, saved_traces_curr, 4, nan);
        saved_PID_traces_all = pad_n_concatenate(saved_PID_traces_all, saved_PID_traces_curr, 4, nan);
        resp_mat_all(:, :, dir_n) = resp_mat;
        
        if suppress_plots == 0
            keyboard
        else
        end
        
        close figure 1
        
    end
    
    del = find(resp_mat_all > 10);
    resp_mat_all(del) = nan;
    del = find(saved_traces_all > 10);
    saved_traces_all(del) = nan;
    
    %plotting averaged response traces
    %paired odor
    color_vecs = [0.65, 0.65, 0.65; 0, 0, 0; 0.65, 0.65, 0.65];
    figure(1)
    
    if set_list_type > 0
        stim_frs_bar = [stim_frs_paired{1}; stim_frs_paired{2}];
        %stim_frs_bar_EL = stim_frs_paired{2};
    elseif set_list_type == 0
        stim_frs_bar = stim_frs_paired{1};
        %stim_frs_bar_EL = stim_frs_paired{2};
    end
    
    for tr_type = 1:2
        mean_trace = squeeze(mean(saved_traces_all(1:(end - 5), 1, tr_type, :), 4, 'omitnan')); 
        se_trace = squeeze(std(saved_traces_all(1:(end - 5), 1, tr_type, :), [], 4, 'omitnan')./sqrt(size(saved_traces_all, 3)));
        col_mult = 0.6.^(tr_type - 1);
        shadedErrorBar([], mean_trace, se_trace, {'Color', color_vecs(tr_type, :)}, 1);
        hold on
    end
    hold off
    ylabel('paired odor responses (dF/F)')
    set_xlabels_time(1, frame_time, 10)
    ax_vals = axis;
    ax_vals(4) = 6;
    ax_vals(3) = -2;
    axis(ax_vals);
    fig_wrapup(1, [])
    if set_list_type == 0
        add_stim_bar(1, stim_frs_bar, paired_color);
    elseif set_list_type == 1
        add_stim_bar(1, stim_frs_bar, [unpaired_color; paired_color]);
    elseif set_list_type == 2
        add_stim_bar(1, stim_frs_bar, [EL_color; paired_color]);
    elseif set_list_type == 3
        add_stim_bar(1, stim_frs_bar, [paired_color; EL_color]);
    else
    end
    
    %unpaired odor
    figure(2)
    for tr_type = 1:2
        mean_trace = squeeze(mean(saved_traces_all(1:(end - 5), 2, tr_type, :), 4, 'omitnan')); 
        se_trace = squeeze(std(saved_traces_all(1:(end - 5), 2, tr_type, :), [], 4, 'omitnan')./sqrt(size(saved_traces_all, 3)));
        col_mult = 0.6.^(tr_type - 1);
        shadedErrorBar([], mean_trace, se_trace, {'Color', color_vecs(tr_type, :)}, 1);
        hold on
    end
    hold off
    ylabel('unpaired odor responses (dF/F)')
    set_xlabels_time(2, frame_time, 10)
    ax_vals = axis;
    ax_vals(4) = 6;
    ax_vals(3) = -2;
    axis(ax_vals);
    fig_wrapup(2, [])
    if set_list_type == 0
        add_stim_bar(2, stim_frs_bar, unpaired_color);
    elseif set_list_type == 1
        add_stim_bar(2, stim_frs_bar, [paired_color; unpaired_color]);
    elseif set_list_type == 2
        add_stim_bar(2, stim_frs_bar, [EL_color; unpaired_color]);
    elseif set_list_type == 3
        add_stim_bar(2, stim_frs_bar, [unpaired_color; EL_color]);
    end

    
    %EL
    figure(3)
    for tr_type = 1:2
        mean_trace = squeeze(mean(saved_traces_all(1:(end - 20), 3, tr_type, :), 4, 'omitnan')); 
        se_trace = squeeze(std(saved_traces_all(1:(end - 20), 3, tr_type, :), [], 4, 'omitnan')./sqrt(size(saved_traces_all, 3)));
        col_mult = 0.6.^(tr_type - 1);
        shadedErrorBar([], mean_trace, se_trace, {'Color', color_vecs(tr_type, :)}, 1);
        hold on
    end
    hold off
    ylabel('EL responses (dF/F)')
    set_xlabels_time(3, frame_time, 10)
    ax_vals = axis;
    ax_vals(4) = 6;
    ax_vals(3) = -2;
    axis(ax_vals);
    fig_wrapup(3, [])
    add_stim_bar(3, stim_frs_bar_EL(1, :), EL_color);
           
    
    %quantification and statistical testing
    
    %comparing pre with post
    %reshaping resp size matrix for plotting function
    resp_mat_small = [];
    if set_list_type <= 1
        n_ods = 3;
    elseif set_list_type > 1
        n_ods = 2;
    else
    end
    
    for odor_n = 1:n_ods
        for resp_type = 1:2
            curr_resps = squeeze(resp_mat_all(resp_type, odor_n, :));
            resp_mat_small = [resp_mat_small, curr_resps];
        end
    end
    
    paired_multiplier = 0.65;
    marker_colors = [paired_color; paired_color.*paired_multiplier; unpaired_color; unpaired_color.*paired_multiplier; EL_color; EL_color.*paired_multiplier];
    marker_colors = marker_colors(1:(n_ods.*2), :);
    line_colors = zeros(size(marker_colors, 1), size(marker_colors, 2)) + 0.7;
    col_pairs = [1, 2; 3, 4; 5, 6];
    col_pairs = col_pairs(1:n_ods, :);
    xlabels = [{'prd pre'}, {'prd post'}, {'unprd pre'}, {'unprd post'}, {'EL pre'}, {'EL post'}];
    xlabels = xlabels(1:(n_ods.*2));
    figure(4)
    fig_h = scattered_dot_plot_ttest(resp_mat_small, 4, 1, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 1, mean_color, 1, 0.05);
    ylabel('response size (dF/F)');
    fig_wrapup(fig_h, []);
    
    
    %comparing paired odor resps with unpaired odor resps
    resp_mat_small_sw = [resp_mat_small(:, 1), resp_mat_small(:, 3), resp_mat_small(:, 2), resp_mat_small(:, 4)];
    marker_colors_sw = [marker_colors(1, :); marker_colors(3, :); marker_colors(2, :); marker_colors(4, :)];
    col_pairs = [1, 2; 3, 4];
    xlabels = [{'prd pre'}, {'unprd pre'}, {'prd post'}, {'unprd post'}];
    figure(5)
    fig_h = scattered_dot_plot_ttest(resp_mat_small_sw, 5, 1, 4, 8, marker_colors_sw, 1, col_pairs, line_colors(1:4, :), xlabels, 1, mean_color, 1, 0.05);
    ylabel('response size (dF/F)');
    fig_wrapup(fig_h, []);
    
    
    
    %comparing post with unlearned
    %reshaping resp size matrix for plotting function
    resp_mat_small = [];
    for odor_n = 1:n_ods
        for resp_type = 2:3
            curr_resps = squeeze(resp_mat_all(resp_type, odor_n, :));
            resp_mat_small = [resp_mat_small, curr_resps];
        end
    end
    
    paired_multiplier = 0.65;
    marker_colors = [ paired_color.*paired_multiplier; paired_color; unpaired_color.*paired_multiplier; unpaired_color; EL_color.*paired_multiplier; EL_color];
    marker_colors = marker_colors(1:(n_ods.*2), :);
    line_colors = marker_colors;
    col_pairs = [1, 2; 3, 4; 5, 6];
    col_pairs = col_pairs(1:n_ods, :);
    xlabels = [{'prd post'}, {'prd unlrn'}, {'unprd post'}, {'unprd unlrn'}, {'EL post'}, {'EL unlrn'}];
    x_labels = xlabels(1:(n_ods.*2));
    figure(6)
    fig_h = scattered_dot_plot_ttest(resp_mat_small, 6, 1, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 1, mean_color, 1, 0.05);
    ylabel('response size (dF/F)');
    fig_wrapup(fig_h, []);
    
    
    %plotting individual traces
    od_type_names = [{'paired '}, {'unpaired '}];
    stim_type_names = [{'pre'}, {'post'}];
    color_vecs = [paired_color; unpaired_color];

    %sorting single traces by total activity in the unpaired, post response
    %dataset
    curr_mat = squeeze(saved_traces_all(:, 2, 2, :))';
    tot_act = sum(curr_mat, 2, 'omitnan');
    
    for od_type = 1:2   %paired-unpaired loop
        for stim_type = 1:2 %pre-post loop
            curr_mat = squeeze(saved_traces_all(:, od_type, stim_type, :))';
            curr_PID_mat = squeeze(saved_PID_traces_all(:, od_type, stim_type, :))';
            
%             %plotting diff traces with mean diff trace
%             fig_h = figure('Name', [od_type_names{od_type}, stim_type_names{stim_type}, ' diff']);
%             curr_mat_diff = diff(curr_mat, 2);
%             plot(curr_mat_diff', 'Color', [0.6, 0.6, 0.6]);
%             hold on
%             %identifying min dif after pulse2 onset
%             [min_amp, min_fr] = min(curr_mat_diff(:, stim_frs_bar(2, 1):(stim_frs_bar(2, 2) + 20)), [], 2);
%             plot((min_fr + stim_frs_bar(2, 1)), min_amp, 'Or');
%             set_xlabels_time(fig_h, frame_time, 10);
%             ylabel('diff. response size');
%             fig_wrapup(fig_h, []);
%             other_type = [1, 2];
%             other_type(other_type == od_type) = [];
%             add_stim_bar(fig_h, stim_frs_bar, [color_vecs(other_type, :); color_vecs(od_type, :)] );
%             
%             %plotting single traces with mean
%             del = find(curr_mat(:, 1) == 1);
%             clust1_mat = curr_mat(del, :);
%             fig_h = figure('Name', [od_type_names{od_type}, stim_type_names{stim_type}]);
%             plot(curr_mat', 'Color', [0.6, 0.6, 0.6]);
%             hold on
%             %plot(mean(curr_mat, 1), 'k', 'lineWidth', 2.5)
%             set_xlabels_time(fig_h, frame_time, 10);
%             ylabel('response size (dF/F)');
%             fig_wrapup(fig_h, []);
%             other_type = [1, 2];
%             other_type(other_type == od_type) = [];
%             add_stim_bar(fig_h, stim_frs_bar, [color_vecs(other_type, :); color_vecs(od_type, :)] );
%             
%             %plotting single PID traces with mean
%             del = find(curr_mat(:, 1) == 1);
%             clust1_mat = curr_mat(del, :);
%             fig_h = figure('Name', [od_type_names{od_type}, stim_type_names{stim_type}, 'PID']);
%             plot(curr_PID_mat', 'Color', [0.6, 0.6, 0.6]);
%             hold on
%             plot(mean(curr_PID_mat, 1), 'k', 'lineWidth', 2.5)
%             set_xlabels_time(fig_h, frame_time, 10);
%             ylabel('response size (dF/F)');
%             fig_wrapup(fig_h, []);
%             other_type = [1, 2];
%             other_type(other_type == od_type) = [];
%             add_stim_bar(fig_h, stim_frs_bar, [color_vecs(other_type, :); color_vecs(od_type, :)] );
%             
                        
            fig_h = figure('Name', [od_type_names{od_type}, stim_type_names{stim_type}, ' single traces']);
            curr_mat = [tot_act, curr_mat];
            curr_mat = sortrows(curr_mat);
            curr_mat(:, 1) = [];
            cascade_plot(fig_h, curr_mat', [0.6, 0.6, 0.6], 1, 0, 1, 3, 3, 1, 1);
            set_xlabels_time(fig_h, frame_time, 10);
            other_type = [1, 2];
            other_type(other_type == od_type) = [];
            %fig_wrapup(fig_h, []);
            add_stim_bar(fig_h, stim_frs_bar, [color_vecs(other_type, :); color_vecs(od_type, :)] );

            
        end
    end
    
    
end