clear all
close all

dataset_list_paths = [...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\MBONG2_PaBaEl_simp_pairing_Berry.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\MBONG2_PaBaEl_handover_pairing_Berry.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\MBONG2_PaBaEl_handover_pairing_Berry_15s_ipi.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\MBONG2_PaBaEl_handover_pairing_Berry_ELsecond.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\MBONG2_PaBaEl_handover_pairing_Berry_longpulse2.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\Berry_handover_MB298B_MBONG4-G1G2_GcaMP6f_starved.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\Berry_handover_13F02_gcaMP6f_starved.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\Berry_handover_GH146_GCaMP6f_fed.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\Berry_handover_SS01240_DPM_starved.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\MBONG2_PaBaEl_handover_pairing_Berry_noLEDctrl.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\Berry_handover_WTG2Ap1_onrig_fedshock.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\Berry_handover_WTG2Ap1_onrig_fedshock_45V.xls'};
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\Berry_handover_WTG2Ap1_onrig_fedshock_30V.xls'};
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\Berry_handover_WTG2Ap1_onrig_fedshock_10V.xls'};
                      
                      ];
            
suppress_plots = 5;
plotting_quant_no_filt = 0;     %1 - only unfiltered traces used for all analysis and plotting - traces included. 0 - filtered traces used for everything.
cell_n = 1;

plot_diff_traces = 0;   %doesn't plot dF/F response traces, but post - pre response traces

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
    
    if isempty(findstr(curr_dir_list_path, 'handover')) == 0
        if isempty(findstr(curr_dir_list_path, 'ELfirst')) == 0
            set_list_type = 2;  %handover stimulus with EL on pulse1 case
        elseif isempty(findstr(curr_dir_list_path, 'ELsecond')) == 0
            set_list_type = 3;  %handover stimulus with EL on pulse2 case
        elseif isempty(findstr(curr_dir_list_path, 'ELfirst')) == 1 && isempty(findstr(curr_dir_list_path, 'ELsecond')) == 1
            set_list_type = 1;  %handover stimulus with no EL pulses
        else
        end
    else 
        set_list_type = 0;  %simple stimulus case
    end
    
    
    saved_traces_all = [];
    saved_PID_traces_all = [];
    ctrast_resps_all = [];
    ctrst_t_win = 1;
    
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
        if dir_n == 1
            ctrst_t_win = round(ctrst_t_win./frame_time);
        else
        end
        
        paired_color = [0,136,55]./256;
        unpaired_color = [166,219,160]./256;
        EL_color = [123,50,148]./256;
        mean_color = [0.8, 0.4, 0.4];
        %mean_color = ([202,0,32]./256);
        
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
        if stim_mat_simple(1, od_olf1_col_n) == 11 & set_list_type ~= 0
            stim_mat(1) = [];
            stim_mat_simple(1, :) = [];
            raw_data_mat(:, :, 1) = [];
            PID_traces(:, 1, :) = [];
        else
        end
        
%         %removing aborted .tif trials
%         l_nan = sum(isnan(raw_data_mat));
%         dummy_trs = find(l_nan(1, 1, :) > 250);
%         raw_data_mat(:, :, dummy_trs) = [];
%         stim_mat(dummy_trs) = [];
%         stim_mat_simple(dummy_trs, :) = [];
%         PID_traces(:, dummy_trs, :) = [];
%         if stim_mat_simple(4, 1) ~= 9
%             dummy_trs = 4;
%             raw_data_mat(:, :, dummy_trs) = [];
%             stim_mat(dummy_trs) = [];
%             stim_mat_simple(dummy_trs, :) = [];
%             PID_traces(:, dummy_trs, :) = [];
%         else
%         end
            
        
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
            
        
        y_ax_lim = [];
        plot_means = 1;
        
        %manually specifying sets of pre and post pairing trial numbers
        if set_list_type <= 1
            tr_lists = [1:3; 7:9; 13:15];   %cases where EL alone trials exist
        elseif set_list_type > 1
            if stim_mat(9).trigger_scan == 0
                tr_lists = [1:2; 6:7; 11:12];   %cases where EL alone trials don't exist
            elseif stim_mat(9).trigger_scan == 1                   
                tr_lists = [1:2; 5:6; 10:11];   %cases where EL alone trials exist
            else
            end
        else
        end
        
        if ctrl_set == 1
            tr_lists = [1:3; 4:6; 7:9];        %no LED stim, so no dummy, pairing trials inserted into stim_mat_simple
        else
        end
       
        
        %looping through pre, post1 and post2 trial types
        resp_mat = zeros(3, 3);
        
        saved_traces_curr = [];
        for tr_type_n = 1:2
            tr_list = tr_lists(tr_type_n, :);
            col_multiplier = 0.6.^(tr_type_n - 1);
            
            %paired odor response
            if set_list_type == 0
                curr_tr = find(stim_mat_simple(tr_list, od_olf1_col_n) == paired_od_n_olf1 & stim_mat_simple(tr_list, led_on_col_n) == 0);
                stim_frs_paired = compute_stim_frs_modular(stim_mat, curr_tr(1), frame_time);
                stim_frs = stim_frs_paired{1};      %integrating over entire pulse1 odor period for simple stim protocols
                stim_frs_pulse1 = stim_frs;
                stim_frs_EL = stim_frs_paired{1};
            elseif set_list_type == 1 || set_list_type == 2     %for handover without EL or handover EL first trials
                curr_tr = find(stim_mat_simple(tr_list, od_olf2_col_n) == paired_od_n_olf2 & stim_mat_simple(tr_list, led_on_col_n) == 0);
                stim_frs_paired = compute_stim_frs_modular(stim_mat, curr_tr(1), frame_time);
                stim_frs_pulse1 = stim_frs_paired{1};
                stim_frs = stim_frs_paired{2};      %integrating over entire pulse2 odor period for simple stim protocols
                stim_frs_EL = stim_frs_paired{1};
            elseif set_list_type == 3       %for handover, EL second trials
                curr_tr = find(stim_mat_simple(tr_list, od_olf1_col_n) == paired_od_n_olf1 & stim_mat_simple(tr_list, led_on_col_n) == 0);
                stim_frs_paired = compute_stim_frs_modular(stim_mat, curr_tr(1), frame_time);
                stim_frs_pulse1 = stim_frs_paired{1};
                stim_frs = stim_frs_paired{2};      %integrating over entire pulse2 odor period for simple stim protocols
                stim_frs_EL = stim_frs_paired{1};
            else
            end
            
            stim_frs(2) = stim_frs(2);  %adding on 2s after odor off to extend integration window
            stim_frs_EL(2) = stim_frs_EL(2); %adding on 2s after odor off to extend integration window
            curr_tr = tr_list(curr_tr);
            
            if length(curr_tr) ~= 1
                disp('more than one trial identified - please disambiguate (paired)')
                disp(['curr_tr = ' num2str(curr_tr)])
                figure(50)
                imagesc(stim_mat_simple, [0, 12]);
                curr_tr = input('Enter correct curr_tr: ');
            else
            end
            size(stim_mat_simple)
        
            curr_trace = squeeze(dff_data_mat_f(:, :, curr_tr));
            resp_vec(1, 1) = mean(mean(curr_trace(stim_frs(1):((stim_frs(1) + integ_win) + round(3./frame_time)) )));
            resp_vec_pulse1(1, 1) = mean(min(curr_trace((stim_frs_pulse1(2) - ctrst_t_win):((stim_frs_pulse1(2))) )));
            ctrst_vec_pulse2(1, 1) = mean(max(curr_trace(stim_frs(1):((stim_frs(1) + ctrst_t_win)) )));
            
            try
                saved_traces_curr(:, 1, tr_type_n) = curr_trace;
            catch
                keyboard
            end
                
            try
                saved_PID_traces_curr(:, 1, tr_type_n) = PID_traces(:, curr_tr);
            catch
                keyboard
            end
                
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
            if length(curr_tr) ~= 1
                disp('more than one trial identified - please disambiguate (unpaired)')
                keyboard
                disp(['curr_tr = ' num2str(curr_tr)])
                figure(50)
                imagesc(stim_mat_simple, [0, 12]);
                curr_tr = input('Enter correct curr_tr: ');
            else
            end
            curr_trace = squeeze(dff_data_mat_f(:, :, curr_tr));
            resp_vec(1, 2) = mean(mean(curr_trace(stim_frs(1):((stim_frs(1) + integ_win) + round(3./frame_time)) )));
            resp_vec_pulse1(1, 2) = mean(min(curr_trace((stim_frs_pulse1(2) - ctrst_t_win):((stim_frs_pulse1(2))) )));
            ctrst_vec_pulse2(1, 2) = mean(max(curr_trace(stim_frs(1):((stim_frs(1) + ctrst_t_win)) )));
            try
                saved_traces_curr(:, 2, tr_type_n) = curr_trace;
            catch
                keyboard
            end
            
            saved_PID_traces_curr(:, 2, tr_type_n) = PID_traces(:, curr_tr);
            subplot(1, 3, 2)
            plot(curr_trace, 'lineWidth', 1.5, 'Color', unpaired_color.*col_multiplier);
            hold on

            %EL odor response
            if set_list_type <= 1
                curr_tr = find(stim_mat_simple(tr_list, od_olf1_col_n) == 11);
                if isempty(curr_tr) == 1
                    curr_tr = find(stim_mat_simple(tr_list, od_olf2_col_n) == 4);
                else
                end
                curr_tr = tr_list(curr_tr);
                
                if length(curr_tr) ~= 1
                    disp('more than one trial identified - please disambiguate (EL)')
                    disp(['curr_tr = ' num2str(curr_tr)])
                    figure(50)
                    imagesc(stim_mat_simple, [0, 12]);
                    curr_tr = input('Enter correct curr_tr: ');
                else
                end
                
                curr_trace = squeeze(dff_data_mat_f(:, :, curr_tr));
                resp_vec(1, 3) = mean(mean(curr_trace(stim_frs_EL(1):((stim_frs_EL(1) + integ_win) + round(3./frame_time)) )));
                
            else
                curr_trace = zeros(size(dff_data_mat_f, 1), 1) + nan;    %since EL alone trials don't exist, padding
                resp_vec(1, 3) = nan;
                
            end
            if set_list_type > 0
                stim_frs_bar_EL = compute_stim_frs_modular(stim_mat, curr_tr(1), frame_time);
                stim_frs_bar_EL = stim_frs_bar_EL{2};
            else
                stim_frs_bar_EL = compute_stim_frs_modular(stim_mat, curr_tr(1), frame_time);
                stim_frs_bar_EL = stim_frs_bar_EL{1};       %for simple stimulus dataset
            end
              
            saved_traces_curr(:, 3, tr_type_n) = curr_trace;
            saved_PID_traces_curr(:, 3, tr_type_n) = PID_traces(:, curr_tr);
            subplot(1, 3, 3)
            plot(curr_trace, 'lineWidth', 1.5, 'Color', EL_color.*col_multiplier);
            hold on
            resp_mat(tr_type_n, :) = resp_vec;
            resp_mat_pulse1(tr_type_n, :) = resp_vec_pulse1;
            ctrst_mat_pulse2(tr_type_n, :) = ctrst_vec_pulse2;
        end
        
        saved_traces_all = pad_n_concatenate(saved_traces_all, saved_traces_curr, 4, nan);
        saved_PID_traces_all = pad_n_concatenate(saved_PID_traces_all, saved_PID_traces_curr, 4, nan);
        clear saved_traces_curr
        clear saved_PID_traces_curr
        resp_mat_all(:, :, dir_n) = resp_mat;
        resp_mat_all_pulse1(:, :, dir_n) = resp_mat_pulse1;
        ctrst_mat_all_pulse2(:, :, dir_n) = ctrst_mat_pulse2;
        
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
        if plot_diff_traces == 1
            diff_traces = saved_traces_all(1:(end - 5), 1, 2, :) - saved_traces_all(1:(end - 5), 1, 1, :);      %post traces - pre traces
            mean_trace = squeeze(mean(diff_traces, 4, 'omitnan'));
            se_trace = squeeze(std(diff_traces, [], 4, 'omitnan')./sqrt(size(saved_traces_all, 3)));
        else
            mean_trace = squeeze(mean(saved_traces_all(1:(end - 5), 1, tr_type, :), 4, 'omitnan')); 
            se_trace = squeeze(std(saved_traces_all(1:(end - 5), 1, tr_type, :), [], 4, 'omitnan')./sqrt(size(saved_traces_all, 3)));
        end
        shadedErrorBar([], mean_trace, se_trace, {'Color', color_vecs(tr_type, :)}, 1);
        if plot_diff_traces == 1
            break
        else
            hold on
        end
    end
    hold off
    if plot_diff_traces == 1
        ylabel('paired odor responses (post - pre)')
    else
        ylabel('paired odor responses (dF/F)')
    end
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
        if plot_diff_traces == 1
            diff_traces = saved_traces_all(1:(end - 5), 2, 2, :) - saved_traces_all(1:(end - 5), 2, 1, :);
            mean_trace = squeeze(mean(diff_traces, 4, 'omitnan')); 
            se_trace = squeeze(std(diff_traces, [], 4, 'omitnan')./sqrt(size(saved_traces_all, 3)));
        else
            mean_trace = squeeze(mean(saved_traces_all(1:(end - 5), 2, tr_type, :), 4, 'omitnan')); 
            se_trace = squeeze(std(saved_traces_all(1:(end - 5), 2, tr_type, :), [], 4, 'omitnan')./sqrt(size(saved_traces_all, 3)));
        end
        
        shadedErrorBar([], mean_trace, se_trace, {'Color', color_vecs(tr_type, :)}, 1);
        if plot_diff_traces == 1
            break
        else
            hold on
        end
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
    
    
    %fitting an exponential decay to mean and individual boosted resps
    if isempty(findstr(curr_dir_list_path, 'longpulse2')) ~= 1
        [del pki] = max(mean_trace(stim_frs(1):(stim_frs(2) + round(5./frame_time)) ));
        pki = pki + stim_frs(1);
        end_pt = min([(stim_frs(2) + round(20./frame_time)), size(mean_trace, 1)]);
        trace_sample = mean_trace(pki:end_pt);
        [x_vec, fit_params] = fit_exp_simple(trace_sample, frame_time);

        modelled_trace = exp(fit_params(2)).*exp(fit_params(1).*x_vec);
        disp(['tau = ',  num2str(-1.*(1./fit_params(1))), ' s']);

        figure(15)
        shadedErrorBar([], mean_trace, se_trace, {'Color', color_vecs(tr_type, :)}, 1);
        hold on
        x_vec = pki:1:(pki + size(modelled_trace, 1) - 1)';
        plot(x_vec, modelled_trace, ':r', 'lineWidth', 2);
        ax_vals = [0, 490, -1, 2];
        axis(ax_vals);
        ylabel('unpaired odor responses (dF/F)')
        set_xlabels_time(15, frame_time, 10)
        fig_wrapup(15, [])
        if set_list_type == 0
            add_stim_bar(15, stim_frs_bar, unpaired_color);
        elseif set_list_type == 1
            add_stim_bar(15, stim_frs_bar, [paired_color; unpaired_color]);
        elseif set_list_type == 2
            add_stim_bar(15, stim_frs_bar, [EL_color; unpaired_color]);
        elseif set_list_type == 3
            add_stim_bar(15, stim_frs_bar, [unpaired_color; EL_color]);
        end
    else
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
    
    %plotting difference traces
    figure(8)
    diff_traces = saved_traces_all(1:(end - 5), 2, 2, :) - saved_traces_all(1:(end - 5), 1, 2, :);      %Unpaired Post - Paired post
    mean_trace = squeeze(mean(diff_traces, 4, 'omitnan'));
    se_trace = squeeze(std(diff_traces, [], 4, 'omitnan')./sqrt(size(saved_traces_all, 3)));
    plot([0, size(mean_trace, 1)], [0, 0], ':', 'lineWidth', 2, 'Color', [0.6, 0.6, 0.6]);
    hold on
    shadedErrorBar([], mean_trace, se_trace, {'Color', color_vecs(tr_type, :)}, 1);
    hold off
    ylabel('unpaired - paired response (dF/F)')
    set_xlabels_time(2, frame_time, 10)
    ax_vals = axis;
    ax_vals(4) = 6;
    ax_vals(3) = -2;
    axis(ax_vals);
    ax_vals = axis();
    
    fig_wrapup(8, [])
    add_stim_bar(8, stim_frs_bar, [0.6, 0.6, 0.6; 0, 0, 0]);
   
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
    
    paired_multiplier = 1;
    marker_colors = [paired_color; paired_color.*paired_multiplier; unpaired_color; unpaired_color.*paired_multiplier; EL_color; EL_color.*paired_multiplier];
    marker_colors = marker_colors(1:(n_ods.*2), :);
    line_colors = zeros(size(marker_colors, 1), size(marker_colors, 2)) + 0.7;
    col_pairs = [1, 2; 3, 4; 5, 6];
    col_pairs = col_pairs(1:n_ods, :);
    xlabels = [{'prd pre'}, {'prd post'}, {'unprd pre'}, {'unprd post'}, {'EL pre'}, {'EL post'}];
    xlabels = xlabels(1:(n_ods.*2));
    figure(4)
    fig_h = scattered_dot_plot_ttest(resp_mat_small, 4, 1, 4, 8, marker_colors, 0, col_pairs, line_colors, xlabels, 1, mean_color, 1, 0.05, 0);
    ylabel('response size (dF/F)');
    fig_wrapup(fig_h, []);
    
    
    %comparing paired odor resps with unpaired odor resps
    resp_mat_small_sw = [resp_mat_small(:, 1), resp_mat_small(:, 3), resp_mat_small(:, 2), resp_mat_small(:, 4)];
    marker_colors_sw = [marker_colors(1, :); marker_colors(3, :); marker_colors(2, :); marker_colors(4, :)];
    col_pairs = [1, 2; 3, 4];
    xlabels = [{'prd pre'}, {'unprd pre'}, {'prd post'}, {'unprd post'}];
    figure(5)
    fig_h = scattered_dot_plot_ttest(resp_mat_small_sw, 5, 1, 4, 8, marker_colors_sw, 0, col_pairs, line_colors(1:4, :), xlabels, 1, mean_color, 1, 0.05, 0);
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
    
    paired_multiplier = 1;
    marker_colors = [ paired_color.*paired_multiplier; paired_color; unpaired_color.*paired_multiplier; unpaired_color; EL_color.*paired_multiplier; EL_color];
    marker_colors = marker_colors(1:(n_ods.*2), :);
    line_colors = marker_colors;
    col_pairs = [1, 2; 3, 4; 5, 6];
    col_pairs = col_pairs(1:n_ods, :);
    xlabels = [{'prd post'}, {'prd unlrn'}, {'unprd post'}, {'unprd unlrn'}, {'EL post'}, {'EL unlrn'}];
    x_labels = xlabels(1:(n_ods.*2));
    figure(6)
    fig_h = scattered_dot_plot_ttest(resp_mat_small, 6, 1, 4, 8, marker_colors, 0, col_pairs, line_colors, xlabels, 1, mean_color, 1, 0.05, 0);
    ylabel('response size (dF/F)');
    fig_wrapup(fig_h, []);
    
    
    
    %Plotting a contrast score (pulse2 - pulse1 response) for transition
    %trials
    contrast_scores = ctrst_mat_all_pulse2 - resp_mat_all_pulse1;
    cscore_mat = [squeeze(contrast_scores(1, 1, :)), squeeze(contrast_scores(1, 2, :)), squeeze(contrast_scores(2, 1, :)), squeeze(contrast_scores(2, 2, :))];
    paired_multiplier = 1;
    marker_colors = [ paired_color; unpaired_color; paired_color.*paired_multiplier;  unpaired_color.*paired_multiplier];
    line_colors = repmat([0.6, 0.6, 0.6], 4, 1);
    col_pairs = [1, 2; 3, 4];
    xlabels = [{'prd pre'}, {'unprd pre'}, {'prd post'}, {'unprd post'}];
    figure(7)
    fig_h = scattered_dot_plot_ttest(cscore_mat, 7, 1, 4, 8, marker_colors, 0, col_pairs, line_colors, xlabels, 1, mean_color, 1, 0.05, 0);
    ylabel('contrast at transition (dF/F)');
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
            if size(curr_mat, 2) > 490
                curr_mat(:, 490:end) = [];
            else
            end
            cascade_plot(fig_h, curr_mat', [0.6, 0.6, 0.6], 1, 0, .3, 3, 3, 0, 1);
            set_xlabels_time(fig_h, frame_time, 10);
            other_type = [1, 2];
            other_type(other_type == od_type) = [];
            %fig_wrapup(fig_h, []);
            
            if set_list_type == 1
                add_stim_bar(fig_h, stim_frs_bar, [color_vecs(other_type, :); color_vecs(od_type, :)] );
            elseif set_list_type == 0
                add_stim_bar(fig_h, stim_frs_bar, [color_vecs(od_type, :)] );
            elseif set_list_type == 2
            add_stim_bar(fig_h, stim_frs_bar, [EL_color; color_vecs(od_type, :)]);
            elseif set_list_type == 3
            add_stim_bar(fig_h, stim_frs_bar, [color_vecs(od_type, :); EL_color]);
            end
            
        end
    end
    
    
end