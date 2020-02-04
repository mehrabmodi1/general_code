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
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_MBONG2_PaBaEl_handover_starved_halfAra_prehabituated_strongUS.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_MBONG2_PaBaEl_handover_starved_halfAra_prehabituated_strongUS_EL_handover.xls'};...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_MBONG2_PaBaEl_handover_starved_halfAra_prehabituated_strongUS_EL_second.xls'};...
                      
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

        %if paired odor is PA, unpaired odor must be BA or vice versa
        unpaired_od_n_olf2 = [1, 3];
        [del, deli] = intersect(unpaired_od_n_olf2, paired_od_n_olf2);
        unpaired_od_n_olf2(deli)= [];
        unpaired_od_name = odor_names2{unpaired_od_n_olf2};
        unpaired_od_n_olf1 = od_name_lookup(odor_names1, unpaired_od_name); 
        
        
        y_ax_lim = [];
        plot_means = 1;
        %plotting and quantifying
        %1. transition trials, paired-unpaired
        [mean_trace_pre, mean_trace_post] = plot_hover_traces(paired_od_n_olf1, unpaired_od_n_olf2, unpaired_od_name, paired_color, unpaired_color, 1,...
            stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, dff_data_mat_f, pairing_tr_n, y_ax_lim, plot_means, suppress_plots);
        
        %2. transition trials, unpaired-paired
        [mean_trace_pre, mean_trace_post] = plot_hover_traces(unpaired_od_n_olf1, paired_od_n_olf2, paired_od_name, unpaired_color, paired_color, 2,...
            stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, dff_data_mat_f, pairing_tr_n, y_ax_lim, plot_means, suppress_plots);
        
        %3. transition trials, EL-paired
        [mean_trace_pre, mean_trace_post] = plot_hover_traces(11, paired_od_n_olf2, paired_od_name, EL_color, paired_color, 3,...
            stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, dff_data_mat_f, pairing_tr_n, y_ax_lim, plot_means, suppress_plots);
        
        %4. transition trials, EL-unpaired
        [mean_trace_pre, mean_trace_post] = plot_hover_traces(11, unpaired_od_n_olf2, unpaired_od_name, EL_color, unpaired_color, 4,...
            stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, dff_data_mat_f, pairing_tr_n, y_ax_lim, plot_means, suppress_plots);
        
        %5. transition trials, paired-EL
        [mean_trace_pre, mean_trace_post] = plot_hover_traces(paired_od_n_olf1, 4, 'Ethyl lactate', paired_color, EL_color, 5,...
            stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, dff_data_mat_f, pairing_tr_n, y_ax_lim, plot_means, suppress_plots);
        
        %6. transition trials, unpaired-EL
        [mean_trace_pre, mean_trace_post] = plot_hover_traces(unpaired_od_n_olf1, 4, 'Ethyl lactate', unpaired_color, EL_color, 6,...
            stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, dff_data_mat_f, pairing_tr_n, y_ax_lim, plot_means, suppress_plots);
        
        
        keyboard
        close all
                
    end
    

end

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