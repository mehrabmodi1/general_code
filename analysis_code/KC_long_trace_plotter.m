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
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_MBONG2_PaBaEl_handover_starved_halfAra_prehabituated_strongUS_EL_second.xls'};...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_c739KC_PaBaEl_handover_prehabituated.xls'};...
                      
                      ];
            
suppress_plots = [];
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
saved_long_traces = 0;
all_sig_frs = [];
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
        curr_dir = manage_base_paths(curr_dir, 3);
        curr_dir = [curr_dir, '\1\'];
       
        tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
        tif_times = tif_times.time_stamps;
        [stim_mat, stim_mat_simple, column_heads, color_vec, good_tr_list, params_orig] = load_params_trains_modular(curr_dir, tif_times);    %reading in trial stimulus parameters after matching time stamps to F traces
        
       
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
        filt_time = .5;            %in s, the time window for boxcar filter for generating filtered traces
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
        od_durs = unique(stim_mat_simple(:, dur_col_ns(2)));
        od_durs(isnan(od_durs)) = [];
        odn_list_olf2 = unique(stim_mat_simple(:, od_col_ns(2)));
        odn_list_olf2(isnan(odn_list_olf2)) = [];
        
        %replacing Nans in stim_mat_simple's od_n_olf2 column with zeroes
        del = isnan(stim_mat_simple(:, od_col_ns(2)));
        del = find(del == 1);
        stim_mat_simple(del, od_col_ns(2)) = 0;
       
        [resp_sizes, sig_trace_mat, sig_trace_mat_old, sig_cell_mat, resp_areaundercurves] = cal_sig_responses_res_modular(dff_data_mat, stim_mat, stim_mat_simple, curr_dir, frame_time, od_col_ns, dur_col_ns);
                
        %KCpop_odor_resp_plotter(dff_data_mat_f, stim_mat_simple, stim_mat, sig_cell_mat, frame_time, 1, [], od_col_ns, dur_col_ns)
        
        %computing mean response trace for each simple odor stimulus presented
        odn_list_olf1 = unique(stim_mat_simple(:, od_col_ns(1)));
        
        for od_n = 1:length(odn_list_olf1)
            od_ni = odn_list_olf1(od_n);
            curr_sig_cells = find(sig_cell_mat(:, od_ni, 2) == 1);
            curr_trs = find(stim_mat_simple(:, od_col_ns(1)) == od_ni & stim_mat_simple(:, od_col_ns(2)) == 0 & stim_mat_simple(:, 22) == 70);
            if isempty(curr_trs) == 1
                continue
            else
            end
            curr_mean_traces = squeeze(mean(dff_data_mat_f(:, curr_sig_cells, curr_trs), 3, 'omitnan'));
            stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
            stim_frs = stim_frs{1};    
            curr_sds = std(curr_mean_traces((stim_frs(1) - round(10./frame_time)):(stim_frs(1) - 2), :), [], 1, 'omitnan');
            curr_means = mean(curr_mean_traces((stim_frs(1) - round(5./frame_time)):(stim_frs(1) - 2), :), 1, 'omitnan');
            curr_threshes = curr_means + 2.33.*curr_sds;
            
            %identifying last point in each response trace that is higher than thresh for each cell.
            for cell_n = 1:length(curr_sig_cells)
                [sig_val_frs, del] =  find(curr_mean_traces(stim_frs(1):end, cell_n) > curr_threshes(cell_n));
                if isempty(sig_val_frs) == 0
                    saved_sig_frs(cell_n, 1) = max(sig_val_frs);
                    if max(sig_val_frs) > 400
                        saved_long_traces = pad_n_concatenate(saved_long_traces, curr_mean_traces(:, cell_n), 2, nan);
                    else
                    end
                else
                    saved_sig_frs(cell_n, 1) = nan;
                end
            end
            all_sig_frs = [all_sig_frs; saved_sig_frs];
        end
                 
        
%         tr_stack = load([curr_dir, 'tr_avg_stack.mat']);
%         tr_stack = tr_stack.ave_stack;
%         playStack(tr_stack, 50, 10)
%         figure(5)
%         plot(saved_sig_frs, '.')
%         keyboard
        close all
                
    end
   
    
    keyboard

end

