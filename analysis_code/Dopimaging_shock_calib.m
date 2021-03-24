clear all
close all

dataset_list_paths = [...
                      
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\Berry_handover_onrig_shockcalib_13F02LexA_DA44.xls'};...
                      
                      ];
            
suppress_plots = 1;
cell_n = 1;

% [del, odor_names1] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
% [del, odor_names2] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList_olf2.xls', 1);
% odor_names2{3} = 'Butyl acetate';

global color_vec;                
a = colormap('bone');
global greymap
greymap = flipud(a);
fly_n = 0;
script_name = mfilename;

elec_color = [0.5, 0.7, 1];

y_ax_traces = 0.8;
y_ax_fit_traces = 0.6;

an_save_path = 'C:\Data\Data\Analysed_data\Analysis_results\PaBaEl_Gamma2\';
force_resave = 1;

n_sec = 2;      %width of moving integration window in s
ROIs_avgd = [3, 6];


integ_win_s = 5;        %width of pulse integration window for response quantification
plotting_quant_no_filt = 0;    %manual toggle to decide whether to use filtered traces for plots, quantification

for list_n = 1:size(dataset_list_paths, 1)
    
    curr_dir_list_path = dataset_list_paths{list_n, 1};
    [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
    n_dirs = size(dir_list, 1);
    dataset_list_name = findstr(curr_dir_list_path, 'list_');
    dataset_list_name = curr_dir_list_path((dataset_list_name + 5):(end - 4));
    
    if isempty(findstr(curr_dir_list_path, 'handover')) == 0
        if isempty(findstr(curr_dir_list_path, 'ELfirst')) == 0
            set_list_type = 2;  %handover stimulus with EL on pulse1 case
        elseif isempty(findstr(curr_dir_list_path, 'ELsecond')) == 0
            set_list_type = 3;  %handover stimulus with EL on pulse1 case
        elseif isempty(findstr(curr_dir_list_path, 'ELfirst')) == 1 && isempty(findstr(curr_dir_list_path, 'ELsecond')) == 1
            set_list_type = 1;  %handover stimulus with no EL pulses
        else
        end
    else 
        set_list_type = 0;  %simple stimulus case
    end
    
    
    saved_traces_all = [];
    saved_PID_traces_all = [];
    saved_resps_all = [];
    saved_odor_resps_all = [];
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
        
        odor_names1 = stim_mat.odourNames;
        odor_names2 = stim_mat.odourNames_olf2;
        PA_odn_olf1 = od_name_lookup(odor_names1, 'Pentyl acetate');
        PA_odn_olf2 = od_name_lookup(odor_names2, 'Pentyl acetate');
        BA_odn_olf1 = od_name_lookup(odor_names1, 'Butyl acetate');
        BA_odn_olf2 = od_name_lookup(odor_names2, 'Butyl acetate');
        
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
%         if cell_n > n_cells
%             cell_n = 1;
%             disp('WARNING! specified cell_n not valid, reassigning cell_n = 1.');
%         else
%         end
        
        %inserting dummy trials (all nans) into raw_data_mat for pairing trials for
        %which no corress .tifs were acquired
        [raw_data_mat, good_tr_list] = match_up_rawmat_matchedtrs(raw_data_mat, stim_mat_simple, tif_n_col_n, good_tr_list);
        
        %if weird EL trial1 exists, getting rid of it
        if stim_mat_simple(1, od_olf1_col_n) == 11 & set_list_type ~= 0
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
          
        if dir_n == 1
            dff_data_mat_f(:, 6, :) = [];
        else
        end
           
            
        
%         if size(dff_data_mat, 2) > 1
%             dff_data_mat = mean(dff_data_mat, 2, 'omitnan');
%             dff_data_mat_f = mean(dff_data_mat_f, 2, 'omitnan');
%         else
%         end
%         dff_data_mat = dff_data_mat(:, cell_n, :);
%         dff_data_mat_f  = dff_data_mat_f(:, cell_n, :);
        del = find(dff_data_mat_f < -1);
        dff_data_mat_f(del) = -1;       %forcing crazy values to sane ones
                
        
        %identifying relevant odor numbers for each olfactometer
        pairing_tr_n = find(stim_mat_simple(:, led_on_col_n) == 1);
        ctrl_set = 0;       %keeping track if this is a no LED ctrl dataset
        if set_list_type == 0
            paired_od_n_olf1 = unique(stim_mat_simple(pairing_tr_n, od_olf1_col_n));    
        elseif set_list_type > 0
            paired_od_n_olf2 = unique(stim_mat_simple(pairing_tr_n, od_olf2_col_n));
            if paired_od_n_olf2 == PA_odn_olf2
                paired_od_n_olf1 = PA_odn_olf1;
            elseif paired_od_n_olf2 == BA_odn_olf2
                paired_od_n_olf1 = BA_odn_olf1;
            elseif isempty(paired_od_n_olf2) == 1       %case for no LED ctrl datasets; paired odor assigned at random
                r_num = rand(1, 1);
                if r_num > 0.5
                    paired_od_n_olf2 = PA_odn_olf2;
                    paired_od_n_olf1 = PA_odn_olf1;
                elseif r_num <= 0.5
                    paired_od_n_olf2 = BA_odn_olf2;
                    paired_od_n_olf1 = BA_odn_olf1;
                else
                end
                ctrl_set = 1;
            end
        else
        end
        
        paired_od_name = odor_names1{paired_od_n_olf1};
      
         
        if paired_od_n_olf1 == PA_odn_olf1
            unpaired_od_n_olf1 = BA_odn_olf1;
            unpaired_od_n_olf2 = BA_odn_olf2;
        elseif paired_od_n_olf1 == BA_odn_olf1
            unpaired_od_n_olf1 = PA_odn_olf1;
            unpaired_od_n_olf2 = PA_odn_olf2;
        else
        end
            
        
        %replacing nans in stim_mat_simple with 0s to make lookup easier
        stim_mat_simple_orig = stim_mat_simple;
        stim_mat_simple(isnan(stim_mat_simple)) = 0;
        
        %identifying shock, no-odor trials
        curr_trs = find(stim_mat_simple(:, 19) == 1 & stim_mat_simple(:, dur_col_ns(1)) < 0.1 & stim_mat_simple(:, dur_col_ns(2)) < 0.1);
        
        curr_traces = mean(dff_data_mat_f(:, ROIs_avgd, curr_trs), 2, 'omitnan');  %averaging across manually specified ROIs
        mean_trace = mean(curr_traces(:, :, 1:3), 3, 'omitnan');   %averaging across repeats
        
        %computing single tr resp sizes
        [stim_frs_train, stim_frs_train_rounded] = compute_LEDelec_stimfrs_train(stim_mat, curr_trs(1), frame_time);
        stim_win = [stim_frs_train_rounded(1, 1), stim_frs_train_rounded(size(stim_frs_train_rounded, 1), size(stim_frs_train_rounded, 2))];
        resps = squeeze(mean(curr_traces(stim_win(1):stim_win(2), :, :), 1));
           
        %logging mean trace for later plotting
        saved_traces_all = pad_n_concatenate(saved_traces_all, squeeze(mean_trace), 2, nan);
        saved_resps_all = [saved_resps_all, squeeze(resps(1:3))];
        
        
        %odor trial analyses
        saved_od_resps = [];
        curr_trs = find(stim_mat_simple(:, od_col_ns(:, 1)) == PA_odn_olf1 & stim_mat_simple(:, od_col_ns(:, 2)) == BA_odn_olf2);   %PA-BA trials
        curr_trs = curr_trs(curr_trs < 4);
        curr_traces = dff_data_mat_f(:, :, curr_trs);  %averaging across manually specified ROIs
        mean_od_resps = mean(curr_traces, 3, 'omitnan');
        saved_od_resps = pad_n_concatenate(saved_od_resps, mean_od_resps, 3, nan);
        
        
        
        curr_trs = find(stim_mat_simple(:, od_col_ns(:, 1)) == BA_odn_olf1 & stim_mat_simple(:, od_col_ns(:, 2)) == PA_odn_olf2);   %BA-PA trials
        curr_trs = curr_trs(curr_trs < 4);
        curr_traces = dff_data_mat_f(:, :, curr_trs);  %averaging across manually specified ROIs
        mean_od_resps = mean(curr_traces, 3, 'omitnan');
        saved_od_resps = pad_n_concatenate(saved_od_resps, mean_od_resps, 3, nan);
        
        curr_trs = find(stim_mat_simple(:, od_col_ns(:, 2)) == 4);   %EL trials
        curr_trs = curr_trs(curr_trs < 4);
        curr_traces = dff_data_mat_f(:, :, curr_trs);  %averaging across manually specified ROIs
        mean_od_resps = mean(curr_traces, 3, 'omitnan');
        saved_od_resps = pad_n_concatenate(saved_od_resps, mean_od_resps, 3, nan);
        
        saved_odor_resps_all = pad_n_concatenate(saved_odor_resps_all, saved_od_resps, 4, nan);
        
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
        figure(3)
        plot(curr_traces)
        ylabel('EL response (dF/F)')
        set_xlabels_time(3, frame_time, 10)
        fig_wrapup(3, []);
        stim_frs_bar = stim_frs{2};
        add_stim_bar(3, stim_frs_bar, [0.6, 0.6, 0.6]);
        
        keyboard
    end
    
    %plotting elec pulse resps
    figure(1)
    plot(saved_traces_all, 'Color', [0.6, 0.6, 0.6])
    hold on
    ylabel('response (dF/F)')
    set_xlabels_time(1, frame_time, 10)
    fig_wrapup(1, []);
    stim_frs_train_rounded(:, 2) = stim_frs_train_rounded(:, 2) + 1;
    add_stim_bar(1, stim_frs_train_rounded, elec_color);
    
    figure(2)
    marker_colors = repmat([0.6, 0.6, 0.6], 3, 1);
    line_colors = zeros(size(marker_colors, 1), size(marker_colors, 2)) + 0.7;
    col_pairs = [];
    xlabels = [{'train1'}, {'train2'}, {'train3'}];
    fig_h = scattered_dot_plot_ttest(saved_resps_all', 4, 1, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 1, mean_color, 1, 0.05);
    ylabel('response (dF/F)');
    fig_wrapup(fig_h, []);
    
end