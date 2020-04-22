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
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_MBONG2_PaBaEl_handover_starved_halfAra_prehabituated_strongUS_EL_second.xls'};...
                      
                      ];
            
suppress_plots = 1;
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

an_save_path = 'C:\Data\Data\Analysed_data\Analysis_results\PaBaEl_Gamma2\';
force_resave = 1;

for list_n = 1:size(dataset_list_paths, 1)
    curr_dir_list_path = dataset_list_paths{list_n, 1};
    [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
    n_dirs = size(dir_list, 1);
    dataset_list_name = findstr(curr_dir_list_path, 'list_');
    dataset_list_name = curr_dir_list_path((dataset_list_name + 5):(end - 4));
    
    save_path = [an_save_path, dataset_list_name, '\'];
    if exist([save_path, 'all_saved_mean_traces_transition_pre.mat']) == 2 && force_resave == 0
        
        all_saved_mean_traces_transition_pre = load([save_path, 'all_saved_mean_traces_transition_pre.mat']); 
        all_saved_mean_traces_transition_pre = all_saved_mean_traces_transition_pre.all_saved_mean_traces_transition_pre;
        all_saved_mean_traces_transition_post = load([save_path, 'all_saved_mean_traces_transition_post.mat']);
        all_saved_mean_traces_transition_post = all_saved_mean_traces_transition_post.all_saved_mean_traces_transition_post;
        all_saved_mean_traces_simple_pre = load([save_path, 'all_saved_mean_traces_simple_pre.mat']);
        all_saved_mean_traces_simple_pre = all_saved_mean_traces_simple_pre.all_saved_mean_traces_simple_pre;
        all_saved_mean_traces_simple_post = load([save_path, 'all_saved_mean_traces_simple_post.mat']);
        all_saved_mean_traces_simple_post = all_saved_mean_traces_simple_post.all_saved_mean_traces_simple_post;
        stim_mat_simple = load([save_path, 'example_stim_mat_simple.mat']);
        stim_mat_simple = stim_mat_simple.stim_mat_simple;
        stim_mat = load([save_path, 'example_stim_mat.mat']);
        stim_mat = stim_mat.stim_mat;
        dur_col_ns = load([save_path, 'example_dur_col_ns.mat']);
        dur_col_ns = dur_col_ns.dur_col_ns;
        frame_time = load([save_path, 'example_frame_time.mat']);
        frame_time = frame_time.frame_time;
        continue
    else
    end

    
    
    %loop to go through all experiment datasets listed in list file
    all_saved_mean_traces_simple_pre = [];
    all_saved_mean_traces_transition_pre = [];
    all_saved_mean_traces_simple_post = [];
    all_saved_mean_traces_transition_post = [];
    
    
    for dir_n = 1:n_dirs
        fly_n = fly_n + 1;
        
        saved_mean_traces_simple_pre = [];
        saved_mean_traces_simple_post = [];
        saved_mean_traces_transition_pre = [];
        saved_mean_traces_transition_post = [];
        
        curr_dir = [dir_list{dir_n, 1}, '\'];
        curr_dir = manage_base_paths(curr_dir, 2);
       
        tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
        tif_times = tif_times.time_stamps;
        [stim_mat, stim_mat_simple, column_heads, color_vec, good_tr_list, params_orig] = load_params_trains_modular(curr_dir, tif_times);    %reading in trial stimulus parameters after matching time stamps to F traces
        paired_color = color_vec(2, :);
        unpaired_color = color_vec(1, :);
        EL_color = color_vec(3, :);
        %mean_color = [0.5, 0.83, 0.98];
        %mean_color = [149, 200, 216]./256;
        mean_color = ([0, 49, 152]./256).*1.5;
        
        stim_mat_simple_nonans = stim_mat_simple;
        stim_mat_simple_nonans(isnan(stim_mat_simple)) = 0;
        
        
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

        %if paired odor is PA, unpaired odor must be BA or vice versa
        unpaired_od_n_olf2 = [1, 3];
        [del, deli] = intersect(unpaired_od_n_olf2, paired_od_n_olf2);
        unpaired_od_n_olf2(deli)= [];
        unpaired_od_name = odor_names2{unpaired_od_n_olf2};
        unpaired_od_n_olf1 = od_name_lookup(odor_names1, unpaired_od_name); 
        
        
        y_ax_lim = [];
        plot_means = 1;
        %plotting and quantifying
                    
        %1. transition trials, unpaired-paired
        [mean_trace_pre, mean_trace_post] = plot_hover_traces(unpaired_od_n_olf1, paired_od_n_olf2, paired_od_name, unpaired_color, paired_color, 2,...
            stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, dff_data_mat_f, pairing_tr_n, y_ax_lim, plot_means, suppress_plots);
        if isempty(mean_trace_pre) == 0
            saved_mean_traces_transition_pre = pad_n_concatenate(saved_mean_traces_transition_pre, mean_trace_pre, 2, nan);
            saved_mean_traces_transition_post = pad_n_concatenate(saved_mean_traces_transition_post, mean_trace_post, 2, nan);
        else
        end
        
        %2. transition trials, paired-unpaired
        [mean_trace_pre, mean_trace_post] = plot_hover_traces(paired_od_n_olf1, unpaired_od_n_olf2, unpaired_od_name, paired_color, unpaired_color, 1,...
            stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, dff_data_mat_f, pairing_tr_n, y_ax_lim, plot_means, suppress_plots);
        if isempty(mean_trace_pre) == 0
            saved_mean_traces_transition_pre = pad_n_concatenate(saved_mean_traces_transition_pre, mean_trace_pre, 2, nan);
            saved_mean_traces_transition_post = pad_n_concatenate(saved_mean_traces_transition_post, mean_trace_post, 2, nan);
        else
        end        
        
        %3. transition trials, EL-paired
        [mean_trace_pre, mean_trace_post] = plot_hover_traces(11, paired_od_n_olf2, paired_od_name, EL_color, paired_color, 3,...
            stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, dff_data_mat_f, pairing_tr_n, y_ax_lim, plot_means, suppress_plots);
        if isempty(mean_trace_pre) == 0
            saved_mean_traces_transition_pre = pad_n_concatenate(saved_mean_traces_transition_pre, mean_trace_pre, 2, nan);
            saved_mean_traces_transition_post = pad_n_concatenate(saved_mean_traces_transition_post, mean_trace_post, 2, nan);
        else
        end
        
        %4. transition trials, EL-unpaired
        [mean_trace_pre, mean_trace_post] = plot_hover_traces(11, unpaired_od_n_olf2, unpaired_od_name, EL_color, unpaired_color, 4,...
            stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, dff_data_mat_f, pairing_tr_n, y_ax_lim, plot_means, suppress_plots);
        if isempty(mean_trace_pre) == 0
            saved_mean_traces_transition_pre = pad_n_concatenate(saved_mean_traces_transition_pre, mean_trace_pre, 2, nan);
            saved_mean_traces_transition_post = pad_n_concatenate(saved_mean_traces_transition_post, mean_trace_post, 2, nan);
        else
        end
        
        %5. transition trials, paired-EL
        [mean_trace_pre, mean_trace_post] = plot_hover_traces(paired_od_n_olf1, 4, 'Ethyl lactate', paired_color, EL_color, 5,...
            stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, dff_data_mat_f, pairing_tr_n, y_ax_lim, plot_means, suppress_plots);
        if isempty(mean_trace_pre) == 0
            saved_mean_traces_transition_pre = pad_n_concatenate(saved_mean_traces_transition_pre, mean_trace_pre, 2, nan);
            saved_mean_traces_transition_post = pad_n_concatenate(saved_mean_traces_transition_post, mean_trace_post, 2, nan);
        else
        end
        
        %6. transition trials, unpaired-EL
        [mean_trace_pre, mean_trace_post] = plot_hover_traces(unpaired_od_n_olf1, 4, 'Ethyl lactate', unpaired_color, EL_color, 6,...
            stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, dff_data_mat_f, pairing_tr_n, y_ax_lim, plot_means, suppress_plots);
        if isempty(mean_trace_pre) == 0
            saved_mean_traces_transition_pre = pad_n_concatenate(saved_mean_traces_transition_pre, mean_trace_pre, 2, nan);
            saved_mean_traces_transition_post = pad_n_concatenate(saved_mean_traces_transition_post, mean_trace_post, 2, nan);
        else
        end
        
        
        %7. simple trials, paired odor
        [mean_trace_pre, mean_trace_post] = plot_simple_traces([], paired_od_n_olf2, paired_od_name, paired_color, 7,...
            stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, dff_data_mat_f, pairing_tr_n, y_ax_lim, plot_means, suppress_plots);
        if isempty(mean_trace_pre) == 0
            saved_mean_traces_simple_pre = pad_n_concatenate(saved_mean_traces_simple_pre, mean_trace_pre, 2, nan);
            saved_mean_traces_simple_post = pad_n_concatenate(saved_mean_traces_simple_post, mean_trace_post, 2, nan);
        else
        end
        
        
        %8. simple trials, unpaired odor
        [mean_trace_pre, mean_trace_post] = plot_simple_traces(unpaired_od_n_olf1, [], unpaired_od_name, unpaired_color, 8,...
            stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, dff_data_mat_f, pairing_tr_n, y_ax_lim, plot_means, suppress_plots);
        if isempty(mean_trace_pre) == 0
            saved_mean_traces_simple_pre = pad_n_concatenate(saved_mean_traces_simple_pre, mean_trace_pre, 2, nan);
            saved_mean_traces_simple_post = pad_n_concatenate(saved_mean_traces_simple_post, mean_trace_post, 2, nan);
        else
        end
        
        %9. simple trials, EL
        [mean_trace_pre, mean_trace_post] = plot_simple_traces(11, [], 'Ethyl lactate', EL_color, 9,...
            stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, dff_data_mat_f, pairing_tr_n, y_ax_lim, plot_means, suppress_plots);
        if isempty(mean_trace_pre) == 0
            saved_mean_traces_simple_pre = pad_n_concatenate(saved_mean_traces_simple_pre, mean_trace_pre, 2, nan);
            saved_mean_traces_simple_post = pad_n_concatenate(saved_mean_traces_simple_post, mean_trace_post, 2, nan);
        else
        end
        
        if suppress_plots == 0
             
        else
        end
        close all
        
        all_saved_mean_traces_transition_pre = pad_n_concatenate(all_saved_mean_traces_transition_pre, saved_mean_traces_transition_pre, 3, nan);
        all_saved_mean_traces_transition_post = pad_n_concatenate(all_saved_mean_traces_transition_post, saved_mean_traces_transition_post, 3, nan);
        
        all_saved_mean_traces_simple_pre = pad_n_concatenate(all_saved_mean_traces_simple_pre, saved_mean_traces_simple_pre, 3, nan);
        all_saved_mean_traces_simple_post = pad_n_concatenate(all_saved_mean_traces_simple_post, saved_mean_traces_simple_post, 3, nan);
                
    end
    

end

%saving response quantifications
save_path = [an_save_path, dataset_list_name, '\'];
if exist([save_path, 'all_saved_mean_traces_transition_pre.mat']) ~= 2 || force_resave == 1
    mkdir(save_path);
    save([save_path, 'all_saved_mean_traces_transition_pre.mat'], 'all_saved_mean_traces_transition_pre');
    save([save_path, 'all_saved_mean_traces_transition_post.mat'], 'all_saved_mean_traces_transition_post');
    save([save_path, 'all_saved_mean_traces_simple_pre.mat'], 'all_saved_mean_traces_simple_pre');
    save([save_path, 'all_saved_mean_traces_simple_post.mat'], 'all_saved_mean_traces_simple_post');
    save([save_path, 'example_stim_mat_simple.mat'], 'stim_mat_simple');
    save([save_path, 'example_stim_mat.mat'], 'stim_mat');
    save([save_path, 'example_dur_col_ns.mat'], 'dur_col_ns');
    save([save_path, 'example_frame_time.mat'], 'frame_time');
else
end

if size(all_saved_mean_traces_simple_pre, 2) == 3
    dataset_type = 2;       %dataset contains EL responses
else
    dataset_type = 1;       %dataset does not contain EL responses
end

%quantification of responses and statistical testing

%quantifying mean responses from different response time-windows
n_sec = 2;      %width of moving integration window in s
transition_trs = find(stim_mat_simple(:, dur_col_ns(1)) == 10 & stim_mat_simple(:, dur_col_ns(2)) == 10);     %list of transition trials
simple_trs = find(stim_mat_simple(:, dur_col_ns(1)) == 10 & stim_mat_simple(:, dur_col_ns(2)) ~= 10 |...
                    stim_mat_simple(:, dur_col_ns(1)) ~= 10 & stim_mat_simple(:, dur_col_ns(2)) == 10);     %list of simple trials

pulse1_on = stim_mat(simple_trs(1)).stimLatency;
pulse1_off = max([stim_mat(simple_trs(1)).duration, stim_mat(simple_trs(1)).duration_olf2])+ pulse1_on;
pulse2_on = stim_mat(transition_trs(1)).rel_stimLatency_olf2 + pulse1_on;
pulse2_off = stim_mat(transition_trs(1)).duration_olf2 + pulse2_on;

stim_frs_trans = [pulse1_on, pulse1_off; pulse2_on, pulse2_off]./frame_time;
stim_frs_simple = [pulse1_on, pulse1_off]./frame_time;
all_integration_wins = [];

%1. integration window = odor period
frames = (1:1:size(all_saved_mean_traces_transition_pre, 1)).*frame_time;
ave_win_simple = (frames > pulse1_on & frames <= pulse1_off);
ave_win_trans = (frames > pulse2_on & frames <= pulse2_off);
whole_od_win_simple_pre = mean(all_saved_mean_traces_simple_pre(ave_win_simple, :, :), 1, 'omitnan');
whole_od_win_simple_post = mean(all_saved_mean_traces_simple_post(ave_win_simple, :, :), 1, 'omitnan');
whole_od_win_trans_pre = mean(all_saved_mean_traces_transition_pre(ave_win_trans, :, :), 1, 'omitnan');
whole_od_win_trans_post = mean(all_saved_mean_traces_transition_post(ave_win_trans, :, :), 1, 'omitnan');
all_integration_wins = [all_integration_wins; ave_win_simple; ave_win_trans];

%2. onset integration window = N s after odor onset
ave_win_simple = (frames > pulse1_on & frames <= (pulse1_on + n_sec) );
ave_win_trans = (frames > pulse2_on & frames <= (pulse2_on + n_sec) );
onset_od_win_simple_pre = mean(all_saved_mean_traces_simple_pre(ave_win_simple, :, :), 1, 'omitnan');
onset_od_win_simple_post = mean(all_saved_mean_traces_simple_post(ave_win_simple, :, :), 1, 'omitnan');
onset_od_win_trans_pre = mean(all_saved_mean_traces_transition_pre(ave_win_trans, :, :), 1, 'omitnan');
onset_od_win_trans_post = mean(all_saved_mean_traces_transition_post(ave_win_trans, :, :), 1, 'omitnan');
all_integration_wins = [all_integration_wins; ave_win_simple; ave_win_trans];

%3. off response integration window = N s after odor off (after pulse2 off for transition trials)
ave_win_simple = (frames > pulse1_off & frames <= (pulse1_off + n_sec) );
ave_win_trans = (frames > pulse2_off & frames <= (pulse2_off + n_sec) );
off_od_win_simple_pre = mean(all_saved_mean_traces_simple_pre(ave_win_simple, :, :), 1, 'omitnan');
off_od_win_simple_post = mean(all_saved_mean_traces_simple_post(ave_win_simple, :, :), 1, 'omitnan');
off_od_win_trans_pre = mean(all_saved_mean_traces_transition_pre(ave_win_trans, :, :), 1, 'omitnan');
off_od_win_trans_post = mean(all_saved_mean_traces_transition_post(ave_win_trans, :, :), 1, 'omitnan');
all_integration_wins = [all_integration_wins; ave_win_simple; ave_win_trans];

%4. extended od period integration window = odor period + N s after odor off
ave_win_simple = (frames > pulse1_on & frames <= (pulse1_off + n_sec) );
ave_win_trans = (frames > pulse2_on & frames <= (pulse2_off + n_sec) );
ext_od_win_simple_pre = mean(all_saved_mean_traces_simple_pre(ave_win_simple, :, :), 1, 'omitnan');
ext_od_win_simple_post = mean(all_saved_mean_traces_simple_post(ave_win_simple, :, :), 1, 'omitnan');
ext_od_win_trans_pre = mean(all_saved_mean_traces_transition_pre(ave_win_trans, :, :), 1, 'omitnan');
ext_od_win_trans_post = mean(all_saved_mean_traces_transition_post(ave_win_trans, :, :), 1, 'omitnan');
all_integration_wins = [all_integration_wins; ave_win_simple; ave_win_trans];

%5 Difference around pulse1-pulse2 transition point, only applies to transition data
ave_win_trans = (frames >= (pulse2_on - n_sec) & frames < (pulse2_on));
pret_od_win_trans_pre = mean(all_saved_mean_traces_transition_pre(ave_win_trans, :, :), 1, 'omitnan');
pret_od_win_trans_post = mean(all_saved_mean_traces_transition_post(ave_win_trans, :, :), 1, 'omitnan');
pret_od_win_trans_pre = onset_od_win_trans_pre - pret_od_win_trans_pre;         %computing difference around pulse1-pulse2 transition
pret_od_win_trans_post = onset_od_win_trans_post - pret_od_win_trans_post;         %computing difference around pulse1-pulse2 transition
all_integration_wins = [all_integration_wins; ave_win_trans];


%Making statistical comparisons

%1. Pre-pairing transition responses are larger than simple trial
%responses
curr_trs = find(stim_mat_simple(:, dur_col_ns(1)) == 10 & stim_mat_simple(:, dur_col_ns(2)) == 10);
trans_od_list = unique(stim_mat_simple(curr_trs, od_col_ns(2)));
fig_n = 1;
figure('Name','Simple v/s Transition odor stim, pre-pairing resps only.')
col_pairs = [1, 2; 3, 4; 5, 6; 7, 8; 9, 10; 11, 12; 13, 14; 15, 16];
if length(trans_od_list) == 2       %case where the paired or unpaired odor are second in transition trials
    onset_resps_paired = [reshape(onset_od_win_simple_pre(:, 1, :), [], 1), reshape(onset_od_win_trans_pre(:, 1, :), [], 1)];
    onset_resps_unpaired = [reshape(onset_od_win_simple_pre(:, 2, :), [], 1), reshape(onset_od_win_trans_pre(:, 2, :), [], 1)];
    off_resps_paired = [reshape(off_od_win_simple_pre(:, 1, :), [], 1), reshape(off_od_win_trans_pre(:, 1, :), [], 1)];
    off_resps_unpaired = [reshape(off_od_win_simple_pre(:, 2, :), [], 1), reshape(off_od_win_trans_pre(:, 2, :), [], 1)];
    whole_resps_paired = [reshape(whole_od_win_simple_pre(:, 1, :), [], 1), reshape(whole_od_win_trans_pre(:, 1, :), [], 1)];
    whole_resps_unpaired = [reshape(whole_od_win_simple_pre(:, 2, :), [], 1), reshape(whole_od_win_trans_pre(:, 2, :), [], 1)];
    ext_resps_paired = [reshape(ext_od_win_simple_pre(:, 1, :), [], 1), reshape(ext_od_win_trans_pre(:, 1, :), [], 1)];
    ext_resps_unpaired = [reshape(ext_od_win_simple_pre(:, 2, :), [], 1), reshape(ext_od_win_trans_pre(:, 2, :), [], 1)];
    plot_mat1 = [onset_resps_paired, onset_resps_unpaired, off_resps_paired, off_resps_unpaired, whole_resps_paired, whole_resps_unpaired, ext_resps_paired, ext_resps_unpaired];
    marker_colors = [paired_color; paired_color; unpaired_color; unpaired_color;...
                        paired_color; paired_color; unpaired_color; unpaired_color;...
                        paired_color; paired_color; unpaired_color; unpaired_color;...
                        paired_color; paired_color; unpaired_color; unpaired_color];
    line_colors = marker_colors;
    
    xlabels = [{'prd_s on'}, {'prd_t on'}, {'unp_s on'}, {'unp_t on'},...
                    {'prd_s off'}, {'prd_t off'}, {'unp_s off'}, {'unp_t off'},...
                    {'prd_s wh'}, {'prd_t wh'}, {'unp_s wh'}, {'unp_t wh'},...
                    {'prd_s ext'}, {'prd_t ext'}, {'unp_s ext'}, {'unp_t ext'}];
    
    
elseif length(trans_od_list) == 1   %case where only EL is the second odor in transition trials
    onset_resps_paired = [reshape(onset_od_win_simple_pre(:, 3, :), [], 1), reshape(onset_od_win_trans_pre(:, 1, :), [], 1)];
    onset_resps_unpaired = [reshape(onset_od_win_simple_pre(:, 3, :), [], 1), reshape(onset_od_win_trans_pre(:, 2, :), [], 1)];
    off_resps_paired = [reshape(off_od_win_simple_pre(:, 3, :), [], 1), reshape(off_od_win_trans_pre(:, 1, :), [], 1)];
    off_resps_unpaired = [reshape(off_od_win_simple_pre(:, 3, :), [], 1), reshape(off_od_win_trans_pre(:, 2, :), [], 1)];
    whole_resps_paired = [reshape(whole_od_win_simple_pre(:, 3, :), [], 1), reshape(whole_od_win_trans_pre(:, 1, :), [], 1)];
    whole_resps_unpaired = [reshape(whole_od_win_simple_pre(:, 3, :), [], 1), reshape(whole_od_win_trans_pre(:, 2, :), [], 1)];
    ext_resps_paired = [reshape(ext_od_win_simple_pre(:, 3, :), [], 1), reshape(ext_od_win_trans_pre(:, 1, :), [], 1)];
    ext_resps_unpaired = [reshape(ext_od_win_simple_pre(:, 3, :), [], 1), reshape(ext_od_win_trans_pre(:, 2, :), [], 1)];
    plot_mat1 = [onset_resps_paired, onset_resps_unpaired, off_resps_paired, off_resps_unpaired, whole_resps_paired, whole_resps_unpaired, ext_resps_paired, ext_resps_unpaired];
    marker_colors = repmat(EL_color, 16, 1);
    line_colors = marker_colors;
    xlabels = [{'EL_s on'}, {'EL_p on'}, {'EL_s on'}, {'EL_u_n_p on'},...
                    {'EL_s off'}, {'EL_p off'}, {'EL_s off'}, {'EL_u_n_p off'},...
                    {'EL_s wh'}, {'EL_p wh'}, {'EL_s wh'}, {'EL_u_n_p wh'},...
                    {'EL_s ext'}, {'EL_p ext'}, {'EL_s ext'}, {'EL_u_n_p ext'},];
else
end
fig_h = scattered_dot_plot_ttest(plot_mat1, fig_n, 1, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 1, mean_color, 1, 0.05);
ylabel('win-averaged response (dF/F)');
fig_wrapup(fig_n, []);


%2. paired_simple_pre vs paired_simple_post, also unpaired and EL as a control
fig_n = fig_n + 1;
figure('Name','Simple pre v/s simple post.')
onset_resps_paired = [reshape(onset_od_win_simple_pre(:, 1, :), [], 1), reshape(onset_od_win_simple_post(:, 1, :), [], 1)];
onset_resps_unpaired = [reshape(onset_od_win_simple_pre(:, 2, :), [], 1), reshape(onset_od_win_simple_post(:, 2, :), [], 1)];
off_resps_paired = [reshape(off_od_win_simple_pre(:, 1, :), [], 1), reshape(off_od_win_simple_post(:, 1, :), [], 1)];
off_resps_unpaired = [reshape(off_od_win_simple_pre(:, 2, :), [], 1), reshape(off_od_win_simple_post(:, 2, :), [], 1)];
whole_resps_paired = [reshape(whole_od_win_simple_pre(:, 1, :), [], 1), reshape(whole_od_win_simple_post(:, 1, :), [], 1)];
whole_resps_unpaired = [reshape(whole_od_win_simple_pre(:, 2, :), [], 1), reshape(whole_od_win_simple_post(:, 2, :), [], 1)];
ext_resps_paired = [reshape(ext_od_win_simple_pre(:, 1, :), [], 1), reshape(ext_od_win_simple_post(:, 1, :), [], 1)];
ext_resps_unpaired = [reshape(ext_od_win_simple_pre(:, 2, :), [], 1), reshape(ext_od_win_simple_post(:, 2, :), [], 1)];

paired_multiplier = 0.7;
if dataset_type == 2
    whole_resps_EL = [reshape(whole_od_win_simple_pre(:, 3, :), [], 1), reshape(whole_od_win_simple_post(:, 3, :), [], 1)];
    EL_colors = [EL_color; EL_color.*paired_multiplier];
    EL_labels = [{'EL_s wh'}, {'EL_s wh'}];
    EL_cols = [17, 18];
else
    whole_resps_EL = [];
    EL_colors = [];
    EL_labels = [];
    EL_cols = [];
end
plot_mat2 = [onset_resps_paired, onset_resps_unpaired, off_resps_paired, off_resps_unpaired, whole_resps_paired, whole_resps_unpaired, ext_resps_paired, ext_resps_unpaired, whole_resps_EL];
plot_mat2z = plot_mat2;
col_pairs = [1, 2; 3, 4; 5, 6; 7, 8; 9, 10; 11, 12; 13, 14; 15, 16; EL_cols];
marker_colors = [paired_color; paired_color.*paired_multiplier; unpaired_color; unpaired_color.*paired_multiplier;...
                    paired_color; paired_color.*paired_multiplier; unpaired_color; unpaired_color.*paired_multiplier;...
                    paired_color; paired_color.*paired_multiplier; unpaired_color; unpaired_color.*paired_multiplier;...
                    paired_color; paired_color.*paired_multiplier; unpaired_color; unpaired_color.*paired_multiplier; EL_colors];
marker_colors2 = marker_colors;
line_colors = marker_colors;

xlabels = [{'prd_s on'}, {'prd_s on'}, {'unp_s on'}, {'unp_s on'},...
                {'prd_s off'}, {'prd_s off'}, {'unp_s off'}, {'unp_s off'},...
                {'prd_s wh'}, {'prd_s wh'}, {'unp_s wh'}, {'unp_s wh'},...
                {'prd_s ext'}, {'prd_s ext'}, {'unp_s ext'}, {'unp_s ext'}, EL_labels];
xlabels2 = xlabels;

fig_h = scattered_dot_plot_ttest(plot_mat2, fig_n, 1, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 1, mean_color, 1, 0.05);
ylabel('win-averaged response (dF/F)');
fig_wrapup(fig_n, []);


%Comparing pre paired resps with pre unpaired resps and the same for post; Simple stimuli.
fig_n = fig_n + 1;
figure('Name','Comparing pre p with pre unp and post p with post unp; Simple');
plot_mat_pre = [whole_resps_paired(:, 1), whole_resps_unpaired(:, 1)];
plot_mat_post = [whole_resps_paired(:, 2), whole_resps_unpaired(:, 2)];
col_pairs = [1, 2];
marker_colors = [paired_color; unpaired_color];
marker_colors_pre = marker_colors.*1.2;
marker_colors_pre(marker_colors_pre > 1) = 1;
mean_color_pre = mean_color.*1.2;
mean_color_pre(mean_color_pre > 1) = 1;
line_colors = repmat([0.65, 0.65, 0.65], 2, 1);

xlabels = [{'paired wh'}, {'unpaired wh'}];

fig_h = scattered_dot_plot_ttest(plot_mat_pre, fig_n, 1, 4, 8, marker_colors_pre, 1, col_pairs, line_colors.*1.2, xlabels, 1, mean_color_pre, 1, 0.05);
hold on
fig_h = scattered_dot_plot_ttest(plot_mat_post, fig_n, 1, 4, 8, marker_colors.*0.6, 1, col_pairs, line_colors.*0.6, xlabels, 1, mean_color.*0.6, 1, 0.05);
hold off
ylabel('win-averaged response (dF/F)');
fig_wrapup(fig_n, []);




%3. paired_trans_pre vs paired_trans_post, also unpaired or if EL post dataset, EL_mean_colortrans_pre vs EL_trans_post
if isempty(strfind(dataset_list_name, 'EL_second')) == 1
    %3. paired_trans_pre vs paired_trans_post
    fig_n = fig_n + 1;
    figure('Name','transition pre v/s transition post.')
    col_pairs = [1, 2; 3, 4; 5, 6; 7, 8; 9, 10; 11, 12; 13, 14; 15, 16];
    onset_resps_paired = [reshape(onset_od_win_trans_pre(:, 1, :), [], 1), reshape(onset_od_win_trans_post(:, 1, :), [], 1)];
    onset_resps_unpaired = [reshape(onset_od_win_trans_pre(:, 2, :), [], 1), reshape(onset_od_win_trans_post(:, 2, :), [], 1)];
    off_resps_paired = [reshape(off_od_win_trans_pre(:, 1, :), [], 1), reshape(off_od_win_trans_post(:, 1, :), [], 1)];
    off_resps_unpaired = [reshape(off_od_win_trans_pre(:, 2, :), [], 1), reshape(off_od_win_trans_post(:, 2, :), [], 1)];
    whole_resps_paired = [reshape(whole_od_win_trans_pre(:, 1, :), [], 1), reshape(whole_od_win_trans_post(:, 1, :), [], 1)];
    whole_resps_unpaired = [reshape(whole_od_win_trans_pre(:, 2, :), [], 1), reshape(whole_od_win_trans_post(:, 2, :), [], 1)];
    ext_resps_paired = [reshape(ext_od_win_trans_pre(:, 1, :), [], 1), reshape(ext_od_win_trans_post(:, 1, :), [], 1)];
    ext_resps_unpaired = [reshape(ext_od_win_trans_pre(:, 2, :), [], 1), reshape(ext_od_win_trans_post(:, 2, :), [], 1)];
    
    plot_mat3 = [onset_resps_paired, onset_resps_unpaired, off_resps_paired, off_resps_unpaired, whole_resps_paired, whole_resps_unpaired, ext_resps_paired, ext_resps_unpaired];
    paired_multiplier = 0.7;
    marker_colors = [paired_color; paired_color.*paired_multiplier; unpaired_color; unpaired_color.*paired_multiplier;...
                        paired_color; paired_color.*paired_multiplier; unpaired_color; unpaired_color.*paired_multiplier;...
                        paired_color; paired_color.*paired_multiplier; unpaired_color; unpaired_color.*paired_multiplier;...
                        paired_color; paired_color.*paired_multiplier; unpaired_color; unpaired_color.*paired_multiplier;...
                        ];
    marker_colors3 = marker_colors;
    line_colors = marker_colors;

    xlabels = [{'prd_t_r on'}, {'prd_t_r on'}, {'unp_t_r on'}, {'unp_t_r on'},...
                    {'prd_t_r off'}, {'prd_t_r off'}, {'unp_t_r off'}, {'unp_t_r off'},...
                    {'prd_t_r wh'}, {'prd_t_r wh'}, {'unp_t_r wh'}, {'unp_t_r wh'},...
                    {'prd_t_r ext'}, {'prd_t_r ext'}, {'unp_t_r ext'}, {'unp_t_r ext'},...
                    {'prd_t_r trdiff'}, {'prd_t_r trdiff'}, {'unp_t_r trdiff'}, {'unp_t_r trdiff'}];
    xlabels3 = xlabels;

    fig_h = scattered_dot_plot_ttest(plot_mat3, fig_n, 1, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 1, mean_color, 1, 0.05);
    ylabel('win-averaged response (dF/F)');
    fig_wrapup(fig_n, []);

else
    %3. EL_trans_pre vs EL_trans_post
    fig_n = fig_n + 1;
    figure('Name','transition pre v/s transition post.')
    onset_resps_paired = [reshape(onset_od_win_trans_pre(:, 1, :), [], 1), reshape(onset_od_win_trans_post(:, 1, :), [], 1)];
    onset_resps_unpaired = [reshape(onset_od_win_trans_pre(:, 2, :), [], 1), reshape(onset_od_win_trans_post(:, 2, :), [], 1)];
    off_resps_paired = [reshape(off_od_win_trans_pre(:, 1, :), [], 1), reshape(off_od_win_trans_post(:, 1, :), [], 1)];
    off_resps_unpaired = [reshape(off_od_win_trans_pre(:, 2, :), [], 1), reshape(off_od_win_trans_post(:, 2, :), [], 1)];
    whole_resps_paired = [reshape(whole_od_win_trans_pre(:, 1, :), [], 1), reshape(whole_od_win_trans_post(:, 1, :), [], 1)];
    whole_resps_unpaired = [reshape(whole_od_win_trans_pre(:, 2, :), [], 1), reshape(whole_od_win_trans_post(:, 2, :), [], 1)];
    ext_resps_paired = [reshape(ext_od_win_trans_pre(:, 1, :), [], 1), reshape(ext_od_win_trans_post(:, 1, :), [], 1)];
    ext_resps_unpaired = [reshape(ext_od_win_trans_pre(:, 2, :), [], 1), reshape(ext_od_win_trans_post(:, 2, :), [], 1)];
    plot_mat3 = [onset_resps_paired, onset_resps_unpaired, off_resps_paired, off_resps_unpaired, whole_resps_paired, whole_resps_unpaired, ext_resps_paired, ext_resps_unpaired];
    paired_multiplier = 0.7;    %for setting paired color darkness
    marker_colors = [EL_color; EL_color.*paired_multiplier; EL_color; EL_color.*paired_multiplier;...
                        EL_color; EL_color.*paired_multiplier; EL_color; EL_color.*paired_multiplier;...
                        EL_color; EL_color.*paired_multiplier; EL_color; EL_color.*paired_multiplier;...
                        EL_color; EL_color.*paired_multiplier; EL_color; EL_color.*paired_multiplier;...
                        EL_color; EL_color.*paired_multiplier; EL_color; EL_color.*paired_multiplier];
    marker_colors3 = marker_colors;
     col_pairs = [1, 2; 3, 4; 5, 6; 7, 8; 9, 10; 11, 12; 13, 14; 15, 16; 17, 18; 19, 20];
    line_colors = marker_colors;

    xlabels = [{'EL_p_r on'}, {'EL_p_r on'}, {'EL_u_n_p on'}, {'EL_u_n_p on'},...
                    {'EL_p_r off'}, {'EL_p_r off'}, {'EL_u_n_p off'}, {'EL_u_n_p off'},...
                    {'EL_p_r wh'}, {'EL_p_r wh'}, {'EL_u_n_p wh'}, {'EL_u_n_p wh'},...
                    {'EL_p_r ext'}, {'EL_p_r ext'}, {'EL_u_n_p ext'}, {'EL_u_n_p ext'},...
                    {'EL_p_r trdiff'}, {'EL_p_r trdiff'}, {'EL_u_n_p trdiff'}, {'EL_u_n_p trdiff'}];
    xlabels3 = xlabels;

    fig_h = scattered_dot_plot_ttest(plot_mat3, fig_n, 1, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 1, mean_color, 1, 0.05);
    ylabel('win-averaged response (dF/F)');
    fig_wrapup(fig_n, []);
end
plot_mat3z = plot_mat3;


%Measuring changes around the odor transition. ie. (pulse2 on response) - (end of pulse1 response)
fig_n = fig_n + 1;
figure('Name','transition pre v/s transition post.')
pret_resps_paired = [reshape(pret_od_win_trans_pre(:, 1, :), [], 1), reshape(ext_od_win_trans_post(:, 1, :), [], 1)];
pret_resps_unpaired = [reshape(pret_od_win_trans_pre(:, 2, :), [], 1), reshape(ext_od_win_trans_post(:, 2, :), [], 1)];

plot_mat = [pret_resps_paired, pret_resps_unpaired];
col_pairs = [1, 2; 3, 4];
paired_multiplier = 0.7;
marker_colors = [paired_color; paired_color.*paired_multiplier; unpaired_color; unpaired_color.*paired_multiplier];
line_colors = marker_colors;

xlabels = [{'prd_t_r'}, {'prd_t_r'}, {'unp_t_r'}, {'unp_t_r'},...
                ];

fig_h = scattered_dot_plot_ttest(plot_mat, fig_n, 1, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 1, mean_color, 1, 0.05);
ylabel('resps around transition (dF/F)');
fig_wrapup(fig_n, []);

%Comparing pre paired resps with pre unpaired resps and the same for post; Transition stimuli.
fig_n = fig_n + 1;
figure('Name','Comparing pre p with pre unp and post p with post unp; Transition');
plot_mat_pre = [onset_resps_paired(:, 1), onset_resps_unpaired(:, 1)];
plot_mat_post = [onset_resps_paired(:, 2), onset_resps_unpaired(:, 2)];
col_pairs = [1, 2];
marker_colors = [paired_color; unpaired_color];
marker_colors_pre = marker_colors.*1.2;
marker_colors_pre(marker_colors_pre > 1) = 1;
mean_color_pre = mean_color.*1.2;
mean_color_pre(mean_color_pre > 1) = 1;
line_colors = repmat([0.65, 0.65, 0.65], 2, 1);
xlabels = [{'paired on'}, {'unpaired on'}];
fig_h = scattered_dot_plot_ttest(plot_mat_pre, fig_n, 1, 4, 8, marker_colors_pre, 1, col_pairs, line_colors.*1.2, xlabels, 1, mean_color_pre, 1, 0.05);
hold on
fig_h = scattered_dot_plot_ttest(plot_mat_post, fig_n, 1, 4, 8, marker_colors.*0.6, 1, col_pairs, line_colors.*0.6, xlabels, 1, mean_color.*0.6, 1, 0.05);
hold off
ylabel('win-averaged response (dF/F)');
fig_wrapup(fig_n, []);
 

%Plotting norm. pre-post differences rather than response sizes: (pre - post)/pre, plotting onset window measures only
%4. simple trial differences
fig_n = fig_n + 1;
figure('Name','norm. Pre-Post responses for simple stimuli, whole window')
plot_mat2i = plot_mat2;
del = find(plot_mat2i(:, 9) < 0.001);
plot_mat2i(del, :) = nan;
del = find(plot_mat2i(:, 11) < 0.001);
plot_mat2i(del, :) = nan;
diff_mat = [(plot_mat2i(:, 9) - plot_mat2i(:, 10))./abs(plot_mat2i(:, 9)), (plot_mat2i(:, 11) - plot_mat2i(:, 12))./abs(plot_mat2i(:, 11))];
paired_multiplier = 0.7;    %for setting paired color darkness
marker_colors = [paired_color; unpaired_color];
col_pairs = [1, 2];
line_colors = repmat([0.6, 0.6, 0.6], 2, 1);
xlabels = [{'paired_s'}, {'unpaired_s'}];
fig_h = scattered_dot_plot_ttest(diff_mat, fig_n, 1, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 1, mean_color, 1, 0.05);
ylabel('Response change (pre-post)/pre');
fig_wrapup(fig_n, []);
 
%5. transition trial differences
fig_n = fig_n + 1;
figure('Name','norm. Pre-Post responses for transition stimuli, on window')
diff_mat2 = [(plot_mat3(:, 1) - plot_mat3(:, 2))./abs(plot_mat3(:, 1)), (plot_mat3(:, 3) - plot_mat2(:, 4))./abs(plot_mat3(:, 3))];
paired_multiplier = 0.7;    %for setting paired color darkness
marker_colors = [paired_color; unpaired_color];
col_pairs = [1, 2];
line_colors = repmat([0.6, 0.6, 0.6], 2, 1);
xlabels = [{'paired_t_r'}, {'unpaired_t_r'}];
fig_h = scattered_dot_plot_ttest(diff_mat2, 5, 1, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 1, mean_color, 1, 0.05);
ylabel('Response change (pre-post)/pre');
fig_wrapup(fig_n, []);
 
%6. testing linear model - subtracting pulse1 response from transition response and re-computing effect of plasticity.
simple_subtracted_paired_pre_traces = all_saved_mean_traces_transition_pre(:, 1, :) - all_saved_mean_traces_simple_pre(:, 2, :);     %subtracting simple unpaired response from transition paired response
simple_subtracted_paired_post_traces = all_saved_mean_traces_transition_post(:, 1, :) - all_saved_mean_traces_simple_post(:, 2, :);     %subtracting simple unpaired response from transition paired response
simple_subtracted_unpaired_pre_traces = all_saved_mean_traces_transition_pre(:, 2, :) - all_saved_mean_traces_simple_pre(:, 1, :);     %subtracting simple unpaired response from transition paired response
simple_subtracted_unpaired_post_traces = all_saved_mean_traces_transition_post(:, 2, :) - all_saved_mean_traces_simple_post(:, 1, :);     %subtracting simple unpaired response from transition paired response

%plotting simple-subtracted transition traces
fig_n = fig_n + 1;
figure('Name','pulse1-subtrtacted paired odor transition responses.')
pairing_multiplier = 0.5;
plot(squeeze(simple_subtracted_paired_pre_traces), 'Color', [0.65, 0.65, 0.65], 'lineWidth', 0.5)
hold on
plot(squeeze(simple_subtracted_paired_post_traces), 'Color', [0, 0, 0], 'lineWidth', 0.5)
plot(squeeze(mean(simple_subtracted_paired_pre_traces, 3, 'omitnan')), 'Color', [0.65, 0.65, 0.65], 'lineWidth', 2)
plot(squeeze(mean(simple_subtracted_paired_post_traces, 3, 'omitnan')), 'Color', [0, 0, 0], 'lineWidth', 2)
hold off
ylabel(['prd od, pulse2 resps (\DeltaF/F)']);
set_xlabels_time(fig_n, frame_time, 10);
fig_wrapup(fig_n, []);
od_color = [0.65, 0.65, 0.65; paired_color];
add_stim_bar(fig_n, stim_frs_trans, od_color);
 
%plotting simple-subtracted transition traces
fig_n = fig_n + 1;
pairing_multiplier = 0.5;
figure('Name','pulse1-subtrtacted unpaired odor transition responses.')
plot(squeeze(simple_subtracted_unpaired_pre_traces), 'Color', [0.65, 0.65, 0.65], 'lineWidth', 0.5)
hold on
plot(squeeze(simple_subtracted_unpaired_post_traces), 'Color', [0, 0, 0], 'lineWidth', 0.5)
plot(squeeze(mean(simple_subtracted_unpaired_pre_traces, 3, 'omitnan')), 'Color', [0.65, 0.65, 0.65], 'lineWidth', 2)
plot(squeeze(mean(simple_subtracted_unpaired_post_traces, 3, 'omitnan')), 'Color', [0, 0, 0], 'lineWidth', 2)
hold off
ylabel(['prd od, pulse2 resps (\DeltaF/F)']);
set_xlabels_time(fig_n, frame_time, 10);
fig_wrapup(fig_n, []);
od_color = [0.65, 0.65, 0.65; unpaired_color];
add_stim_bar(fig_n, stim_frs_trans, od_color);

%computing responses
fig_n = fig_n + 1;
figure('Name','Transition responses after subtracting simple, pulse1 response.')
ave_win_trans = (frames > pulse2_on & frames <= (pulse2_on + n_sec) );
sim_subtr_paired_pre = squeeze(mean(simple_subtracted_paired_pre_traces(ave_win_trans, :, :), 1, 'omitnan'));
sim_subtr_paired_post = squeeze(mean(simple_subtracted_paired_post_traces(ave_win_trans, :, :), 1, 'omitnan'));
sim_subtr_unpaired_pre = squeeze(mean(simple_subtracted_unpaired_pre_traces(ave_win_trans, :, :), 1, 'omitnan'));
sim_subtr_unpaired_post = squeeze(mean(simple_subtracted_unpaired_post_traces(ave_win_trans, :, :), 1, 'omitnan'));

plot_mat = [sim_subtr_paired_pre, sim_subtr_paired_post, sim_subtr_unpaired_pre, sim_subtr_unpaired_post];

paired_multiplier = 0.7;
marker_colors = [paired_color; paired_color.*paired_multiplier; unpaired_color; unpaired_color.*paired_multiplier];
                   
line_colors = marker_colors;
col_pairs = [1, 2; 3, 4];
xlabels = [{'prd_s_u_b pre'}, {'prd_s_u_b post'}, {'unprd_s_u_b pre'}, {'unprd_s_u_b post'}];

fig_h = scattered_dot_plot_ttest(plot_mat, fig_n, 1, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 1, mean_color, 1, 0.05);
ylabel('onset win-averaged response (dF/F)');
fig_wrapup(fig_n, []);


%7 Plotting fly-averaged traces mean+-SE
%simple traces
mean_paired_pre = mean(squeeze(all_saved_mean_traces_simple_pre(:, 1, :)), 2, 'omitnan');
se_paired_pre = std(squeeze(all_saved_mean_traces_simple_pre(:, 1, :)), [], 2, 'omitnan')./sqrt(size(all_saved_mean_traces_simple_pre, 3));
mean_unpaired_pre = mean(squeeze(all_saved_mean_traces_simple_pre(:, 2, :)), 2, 'omitnan');
se_unpaired_pre = std(squeeze(all_saved_mean_traces_simple_pre(:, 2, :)), [], 2, 'omitnan')./sqrt(size(all_saved_mean_traces_simple_pre, 3));
max_y = max((mean_paired_pre + se_paired_pre));
max_y = max([max_y; max((mean_unpaired_pre + se_unpaired_pre))]);

mean_paired_post = mean(squeeze(all_saved_mean_traces_simple_post(:, 1, :)), 2, 'omitnan');
se_paired_post = std(squeeze(all_saved_mean_traces_simple_post(:, 1, :)), [], 2, 'omitnan')./sqrt(size(all_saved_mean_traces_simple_post, 3));
mean_unpaired_post = mean(squeeze(all_saved_mean_traces_simple_post(:, 2, :)), 2, 'omitnan');
se_unpaired_post = std(squeeze(all_saved_mean_traces_simple_post(:, 2, :)), [], 2, 'omitnan')./sqrt(size(all_saved_mean_traces_simple_post, 3));
max_y = max([max_y; max((mean_paired_post + se_paired_post))]);
max_y = max([max_y; max((mean_unpaired_post + se_unpaired_post))]);

if dataset_type == 2
    mean_EL_pre = mean(squeeze(all_saved_mean_traces_simple_pre(:, 3, :)), 2, 'omitnan');
    se_EL_pre = std(squeeze(all_saved_mean_traces_simple_pre(:, 3, :)), [], 2, 'omitnan')./sqrt(size(all_saved_mean_traces_simple_pre, 3));
    max_y = max([max_y; max((mean_EL_pre + se_EL_pre))]);
    mean_EL_post = mean(squeeze(all_saved_mean_traces_simple_post(:, 3, :)), 2, 'omitnan');
    se_EL_post = std(squeeze(all_saved_mean_traces_simple_post(:, 3, :)), [], 2, 'omitnan')./sqrt(size(all_saved_mean_traces_simple_post, 3));
    max_y = max([max_y; max((mean_EL_post + se_EL_post))]);
    
else
end


%transition traces
mean_paired_pre_t = mean(squeeze(all_saved_mean_traces_transition_pre(:, 1, :)), 2, 'omitnan');
se_paired_pre_t = std(squeeze(all_saved_mean_traces_transition_pre(:, 1, :)), [], 2, 'omitnan')./sqrt(size(all_saved_mean_traces_transition_pre, 3));
mean_unpaired_pre_t = mean(squeeze(all_saved_mean_traces_transition_pre(:, 2, :)), 2, 'omitnan');
se_unpaired_pre_t = std(squeeze(all_saved_mean_traces_transition_pre(:, 2, :)), [], 2, 'omitnan')./sqrt(size(all_saved_mean_traces_transition_pre, 3));
max_y = max([max_y; max((mean_paired_pre_t + se_paired_pre_t))]);
max_y = max([max_y; max((mean_unpaired_pre_t + se_unpaired_pre_t))]);

mean_paired_post_t = mean(squeeze(all_saved_mean_traces_transition_post(:, 1, :)), 2, 'omitnan');
se_paired_post_t = std(squeeze(all_saved_mean_traces_transition_post(:, 1, :)), [], 2, 'omitnan')./sqrt(size(all_saved_mean_traces_transition_post, 3));
mean_unpaired_post_t = mean(squeeze(all_saved_mean_traces_transition_post(:, 2, :)), 2, 'omitnan');
se_unpaired_post_t = std(squeeze(all_saved_mean_traces_transition_post(:, 2, :)), [], 2, 'omitnan')./sqrt(size(all_saved_mean_traces_transition_post, 3));
max_y = max([max_y; max((mean_paired_post_t + se_paired_post_t))]);
max_y = max([max_y; max((mean_unpaired_post_t + se_unpaired_post_t))]);


fig_n = fig_n + 1;
figure('Name','Simple responses, paired od.')
shadedErrorBar([], mean_paired_pre, se_paired_pre, {'Color', [0.65, 0.65, 0.65], 'lineWidth', 2}, 1);
hold on
shadedErrorBar([], mean_paired_post, se_paired_post, {'Color', [0, 0, 0], 'lineWidth', 2}, 1);
ax_vals = axis;
ax_vals(4) = max_y.*1.2;
axis(ax_vals);
ylabel('dF/F response');
set_xlabels_time(fig_n, frame_time, 10);
fig_wrapup(fig_n, []);
od_color = [paired_color];
add_stim_bar(fig_n, stim_frs_simple, od_color);


fig_n = fig_n + 1;
figure('Name','Simple responses, unpaired od.')
shadedErrorBar([], mean_unpaired_pre, se_unpaired_pre, {'Color', [0.65, 0.65, 0.65], 'lineWidth', 2}, 1);
hold on
shadedErrorBar([], mean_unpaired_post, se_unpaired_post, {'Color', [0, 0, 0], 'lineWidth', 2}, 1);
ax_vals = axis;
ax_vals(4) = max_y.*1.2;
axis(ax_vals);
ylabel('dF/F response');
set_xlabels_time(fig_n, frame_time, 10);
fig_wrapup(fig_n, []);
od_color = [unpaired_color];
add_stim_bar(fig_n, stim_frs_simple, od_color);


fig_n = fig_n + 1;
figure('Name','Simple responses, EL.')
if dataset_type == 2
    shadedErrorBar([], mean_EL_pre, se_EL_pre, {'Color', [0.65, 0.65, 0.65], 'lineWidth', 2}, 1);
    hold on
    shadedErrorBar([], mean_EL_post, se_EL_post, {'Color', [0, 0, 0], 'lineWidth', 2}, 1);
    ax_vals = axis;
    ax_vals(4) = max_y.*1.2;
    axis(ax_vals);
    ylabel('dF/F response');
    set_xlabels_time(fig_n, frame_time, 10);
    fig_wrapup(fig_n, []);
    od_color = [EL_color];
    add_stim_bar(fig_n, stim_frs_simple, od_color);
else
end

fig_n = fig_n + 1;
figure('Name','Transition responses, paired od.')
shadedErrorBar([], mean_paired_pre_t, se_paired_pre_t, {'Color', [0.65, 0.65, 0.65], 'lineWidth', 2}, 1);
hold on
shadedErrorBar([], mean_paired_post_t, se_paired_post_t, {'Color', [0, 0, 0], 'lineWidth', 2}, 1);
ax_vals = axis;
ax_vals(4) = max_y.*1.2;
axis(ax_vals);
ylabel('dF/F response');
set_xlabels_time(fig_n, frame_time, 10);
fig_wrapup(fig_n, []);
od_color = [[0.65, 0.65, 0.65]; paired_color];
add_stim_bar(fig_n, stim_frs_trans, od_color);


fig_n = fig_n + 1;
figure('Name','Transition responses, unpaired od.')
shadedErrorBar([], mean_unpaired_pre_t, se_unpaired_pre_t, {'Color', [0.65, 0.65, 0.65], 'lineWidth', 2}, 1);
hold on
shadedErrorBar([], mean_unpaired_post_t, se_unpaired_post_t, {'Color', [0, 0, 0], 'lineWidth', 2}, 1);
ax_vals = axis;
ax_vals(4) = max_y.*1.2;
axis(ax_vals);
ylabel('dF/F response');
set_xlabels_time(fig_n, frame_time, 10);
fig_wrapup(fig_n, []);
od_color = [[0.65, 0.65, 0.65]; unpaired_color];
add_stim_bar(fig_n, stim_frs_trans, od_color);


%PICK UP THREAD HERE
%figure out discrepancy between these figs and figure 4. Also, separate
%last two columns from figure 4 into a separate plot.

%8 Plotting simple bar plots to match Berry et al.
%simple bar plot
fig_n = fig_n + 1;
figure('Name','Simple pre v/s simple post. Whole int window')
col_list = 9:1:12;
if dataset_type == 2
    col_list = [col_list, 17, 18];
else
    
end

[p_vec_out, effect_sizes_out, min_dif_out] = stat_testing(plot_mat2(:, col_list), 1, 0, 0.05);
stats = [p_vec_out, effect_sizes_out, min_dif_out];
mean_vec = mean(plot_mat2(:, col_list), 1);
se_vec = std(plot_mat2(:, col_list), [], 1)./sqrt(size(plot_mat2(:, 9:12), 1));
plot_vec = 1:1:length(mean_vec);
b = bar(plot_vec, mean_vec);
b.FaceColor = 'flat';

%inserting text labels for p values etc
x_pos = 1.5:2:3.5;
y_pos = [max(mean_vec(1:2)), max(mean_vec(3:4))].*1.2;
add_plot_labels(fig_n, stats, x_pos, y_pos);


for bar_n = 1:length(col_list)
     colorn = col_list(bar_n);
    b.CData(bar_n, :) = marker_colors2(colorn, :);
end
hold on
errorbar(plot_vec, mean_vec, se_vec, 'LineStyle', 'none', 'lineWidth', 1.5, 'Color', 'k');
ylabel('win-averaged response (dF/F)');
ax = gca;
ax.XTick = plot_vec;
ax.XTickLabels = xlabels2(col_list);
fig_wrapup(fig_n, []);
hold off


%transition bar plot
fig_n = fig_n + 1;
figure('Name','Trans pre v/s trans post. On int window')

col_list = 1:1:4;
if dataset_type == 2
    col_list = [col_list, 17, 18];
else
    
end

[p_vec_out, effect_sizes_out, min_dif_out] = stat_testing(plot_mat3(:, col_list), 1, 0, 0.05);
stats = [p_vec_out, effect_sizes_out, min_dif_out];
mean_vec = mean(plot_mat3(:, 1:4), 1);
se_vec = std(plot_mat3(:, 1:4), [], 1)./sqrt(size(plot_mat3(:, 1:4), 1));
plot_vec = 1:1:length(mean_vec);
b = bar(plot_vec, mean_vec);
%inserting text labels for p values etc
x_pos = 1.5:2:3.5;
y_pos = [max(mean_vec(1:2)), max(mean_vec(3:4))].*1.2;
add_plot_labels(fig_n, stats, x_pos, y_pos);
b.FaceColor = 'flat';

for bar_n = 1:4
    b.CData(bar_n, :) = marker_colors3(bar_n, :);
end
hold on
errorbar(plot_vec, mean_vec, se_vec, 'LineStyle', 'none', 'lineWidth', 1.5, 'Color', 'k');
ylabel('win-averaged response (dF/F)');
ax = gca;
ax.XTick = plot_vec;
ax.XTickLabels = xlabels3(1:4);
fig_wrapup(fig_n, []);
hold off


%---------------------------------
%worker functions
function [mean_trace_pre, mean_trace_post] = plot_hover_traces(olf1_odn, olf2_odn, od2_name, od1_color, od2_color, fign,...
    stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, dff_data_mat_f, pairing_tr, y_ax_lim, plot_means, suppress_plots)
frame_time = 0.099;

%identifying current trials
curr_trs = find(stim_mat_simple(:, od_col_ns(1)) == olf1_odn & stim_mat_simple(:, dur_col_ns(1)) == 10 ...
    & stim_mat_simple(:, od_col_ns(2)) == olf2_odn & stim_mat_simple(:, dur_col_ns(2)) == 10);

%case where current od1-od2 transition isn't part of the dataset
if isempty(curr_trs) == 1
    mean_trace_pre = [];
    mean_trace_post = [];
    return
else
end

curr_trs_pre = curr_trs(curr_trs < pairing_tr);
curr_trs_post = curr_trs(curr_trs > pairing_tr);

curr_traces_pre = squeeze(dff_data_mat_f(:, :, curr_trs_pre));
curr_traces_post = squeeze(dff_data_mat_f(:, :, curr_trs_post));
mean_trace_pre = mean(curr_traces_pre, 2, 'omitnan');
mean_trace_post = mean(curr_traces_post, 2, 'omitnan');

if suppress_plots == 0
    stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
    stim_frs = [stim_frs{1}; stim_frs{2}]; 
    figure(fign)

    if plot_means == 1
        plot(mean_trace_pre, 'lineWidth', 2.5, 'Color', [0.65, 0.65, 0.65]);
        hold on
        plot(mean_trace_post, 'lineWidth', 2.5, 'Color', [0, 0, 0]);
        plot(curr_traces_pre, 'lineWidth', 0.5, 'Color', [0.65, 0.65, 0.65]);
        plot(curr_traces_post, 'lineWidth', 0.5, 'Color', [0, 0, 0]);
    elseif plot_means == 0
        plot(curr_traces_pre, 'lineWidth', 0.5, 'Color', [0.65, 0.65, 0.65]);
        hold on
        plot(curr_traces_post, 'lineWidth', 0.5, 'Color', [0, 0, 0]);
    else
    end
    hold off
    ylabel([od2_name, ' responses (\DeltaF/F)']);
    
    if isempty(y_ax_lim) == 0
        ax_vals = axis;
        ax_vals(4) = y_ax_lim;
        axis(ax_vals);
    else
    end
    
    set_xlabels_time(fign, frame_time, 10);
    fig_wrapup(fign, []);
    add_stim_bar(fign, stim_frs, [od1_color; od2_color]);

elseif suppress_plots == 1
else
end

end



function [mean_trace_pre, mean_trace_post] = plot_simple_traces(olf1_odn, olf2_odn, od_name, od_color, fign,...
    stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, dff_data_mat_f, pairing_tr, y_ax_lim, plot_means, suppress_plots)
frame_time = 0.099;

if isempty(olf1_odn) == 1
    olf_n = 2;
    olf_n_other = 1;
    od_n = olf2_odn;
elseif isempty(olf2_odn) == 1
    olf_n = 1;
    olf_n_other = 2;
    od_n = olf1_odn;
else
end

stim_mat_simple_nonans = stim_mat_simple;
stim_mat_simple_nonans(isnan(stim_mat_simple_nonans)) = 0;

%identifying current trials
curr_trs = find(stim_mat_simple(:, od_col_ns(olf_n)) == od_n & stim_mat_simple(:, dur_col_ns(olf_n)) == 10 ...
    & stim_mat_simple_nonans(:, dur_col_ns(olf_n_other)) < 1);

%case where current od1-od2 transition isn't part of the dataset
if isempty(curr_trs) == 1
    mean_trace_pre = [];
    mean_trace_post = [];
    return
else
end

curr_trs_pre = curr_trs(curr_trs < pairing_tr);
curr_trs_post = curr_trs(curr_trs > pairing_tr);

curr_traces_pre = squeeze(dff_data_mat_f(:, :, curr_trs_pre));
curr_traces_post = squeeze(dff_data_mat_f(:, :, curr_trs_post));
mean_trace_pre = mean(curr_traces_pre, 2, 'omitnan');
mean_trace_post = mean(curr_traces_post, 2, 'omitnan');

if suppress_plots == 0
    stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
    stim_frs = stim_frs{olf_n}; 
    figure(fign)

    if plot_means == 1
        plot(mean_trace_pre, 'lineWidth', 2.5, 'Color', [0.65, 0.65, 0.65]);
        hold on
        plot(mean_trace_post, 'lineWidth', 2.5, 'Color', [0, 0, 0]);
        plot(curr_traces_pre, 'lineWidth', 0.5, 'Color', [0.65, 0.65, 0.65]);
        plot(curr_traces_post, 'lineWidth', 0.5, 'Color', [0, 0, 0]);
    elseif plot_means == 0
        plot(curr_traces_pre, 'lineWidth', 0.5, 'Color', [0.65, 0.65, 0.65]);
        hold on
        plot(curr_traces_post, 'lineWidth', 0.5, 'Color', [0, 0, 0]);
    else
    end
    hold off
    ylabel([od_name, ' responses (\DeltaF/F)']);
    
    if isempty(y_ax_lim) == 0
        ax_vals = axis;
        ax_vals(4) = y_ax_lim;
        axis(ax_vals);
    else
    end
    
    set_xlabels_time(fign, frame_time, 10);
    fig_wrapup(fign, []);
    add_stim_bar(fign, stim_frs, od_color);

elseif suppress_plots == 1
else
end

end


%function to add labels to plots for p values etc.
function [] = add_plot_labels(fig_n, stats, x_pos, y_pos)
    figure(fig_n);
    n_pairs = size(stats, 1);
    for pair_n = 1:n_pairs
        curr_p_val = round(stats(pair_n, 1), 3);
        curr_d = round(stats(pair_n, 2), 3);
        curr_dif = round(stats(pair_n, 3), 3);
        curr_ac_dif = round(stats(pair_n, 4), 3);
        min_n = stats(pair_n, 5);
        ac_n = stats(pair_n, 6);
        
         if curr_p_val ~= 0
                p_label = ['p = ', num2str(curr_p_val)];
            elseif curr_p_val == 0
                p_label = ['p < ', '0.001'];
            else
            end
            if curr_d ~= 0
                d_label = ['Cohen''s d = ', num2str(curr_d)];
            elseif curr_d == 0
                d_label = ['Cohen''s d  < ', '0.001'];
            else
            end
            n_label = ['min n = ' num2str(min_n), ', n = ' num2str(ac_n)];
            if curr_dif ~= 0
                dif_label = ['min. dif = ', num2str(curr_dif), ', dif = ' num2str(curr_ac_dif)];
            elseif curr_dif == 0
                dif_label = ['min. dif  < ', '0.001, dif = ', num2str(curr_ac_dif)];
            elseif isempty(curr_dif) == 1
                dif_label = [];
                n_label = [];
            else
            end
            x_posi = x_pos(pair_n);
            y_posi = y_pos(pair_n);
            lin_spac = max(y_pos).*0.055;
            text(x_posi, y_posi, p_label);
            text(x_posi, (y_posi + lin_spac), d_label);
%             text(x_posi, (y_posi + 2.*lin_spac), dif_label);
%             text(x_posi, (y_posi + 3.*lin_spac), n_label);
        
    end
   
end
