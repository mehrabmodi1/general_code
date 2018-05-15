clear all
close all

[del, odor_names] = xlsread('C:\Data\Code\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
dir_list_path = 'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_odpulse_freq_series.xls';

[del, dir_list] = xlsread(dir_list_path, 1);        %list of Suite2P results directories
n_dirs = size(dir_list, 1);

a = colormap('bone');
global greymap
greymap = flipud(a);


%loop to go through all experiment datasets listed in list file
for dir_n = 1:n_dirs
    curr_dir = [dir_list{dir_n, 1}, '\'];
    tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
    tif_times = tif_times.time_stamps;
    [stim_mat, stim_mat_simple, column_heads, color_vec] = load_params_trains(curr_dir, tif_times);    %reading in trial stimulus parameters after matching time stamps to F traces
    
    %Reading in experimental parameters
    odor_list = unique(stim_mat_simple(:, 2) );
    n_odors = length(odor_list);
    odor_dur_list = unique(stim_mat_simple(:, 3) );
    n_od_durs = length(odor_dur_list);
    n_trains = max(stim_mat_simple(:, 11));
    
    cd(curr_dir);
    tif_name = dir('*.tif');
    stack_obj = ScanImageTiffReader([curr_dir, tif_name(1).name]);
    [frame_time, zoom, n_chans] = SI_tif_info(stack_obj);
    
    %loading extracted raw fluorescence data matrices written by raw_dff_extractor
    raw_data_mat = load([curr_dir 'extracted_raw_data_mat.mat']);
    raw_data_mat = raw_data_mat.raw_data_mat;           %raw F traces extracted from ROIs
    raw_data_mat_orig = raw_data_mat;
    raw_data_mat = raw_data_mat(:, :, stim_mat_simple(:, 1));       %making sure only time-stamp matched trials are used for further analysis

    %calculating dF/F traces from raw data
    filt_time = 0.5;            %in s, the time window for boxcar filter for generating filtered traces
    [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, curr_dir);

    %identifying significantly responsive cells
    [resp_areas, sig_trace_mat, sig_trace_mat_old, sig_cell_mat] = cal_sig_responses_res(dff_data_mat, stim_mat, stim_mat_simple, curr_dir, frame_time);

    %Running data quality control checks
    sig_cell_mat_old = sig_cell_mat;
    [sig_cell_mat, all_bad_trs] = cell_data_quality_control(dff_data_mat_f, stim_mat, stim_mat_simple, sig_cell_mat, 1, frame_time);
    dff_data_mat(:, :, all_bad_trs) = nan;
    %disp([num2str((length(all_bad_trs)./size(dff_data_mat, 3)).*100) ' percent of trials were auto-identified as bad and removed.']);
    sig_cells = find(sum(sum(sig_cell_mat, 3), 2) > 0);         %list of all cells significant for any odor for any duration
    
    %% Analysing response data

    %plotting ave resp mats of sig cells for different durs
    for odor_n = 1:n_odors
        odor_ni = odor_list(odor_n);        %actual odor number rather than list count
        
        for dur_n = 1:n_od_durs
            curr_dur = odor_dur_list(dur_n);
            curr_trs = find(stim_mat_simple(:, 2) == odor_ni & stim_mat_simple(:, 3) == curr_dur);
            curr_sig_cells = find(sig_cell_mat(:, odor_ni, dur_n) == 1);
            curr_insig_cells = find(sig_cell_mat(:, odor_ni, dur_n) == 0);
            ave_mat = mean(dff_data_mat_f(:, curr_sig_cells, curr_trs), 3, 'omitnan');
            ave_mat_insig = mean(dff_data_mat_f(:, curr_insig_cells, curr_trs), 3, 'omitnan');
            stim_frs = compute_stim_frs(stim_mat, curr_trs(1), frame_time);
            
            
            %plotting
            figure(1)
            imagesc(ave_mat')
            colormap(greymap)
            add_stim_bar(1, stim_frs, color_vec(odor_n, :));
            title('sig responders')
            
            figure(2)
            imagesc(ave_mat_insig')
            colormap(greymap)
            add_stim_bar(2, stim_frs, color_vec(odor_n, :));
            title('insig responders')
            keyboard
            close figure 1
            close figure 2
            
            
        end
    end
    
    
    
    %looking for exclusive responders ie. cell-od pairs responsive to only one duration
    
%     for odor_n = 1:n_odors
%         odor_ni = odor_list(odor_n);        %actual odor number rather than list count
%         
%         for dur_n = 1:n_od_durs
%             
% 
%         end
%     end
    
    keyboard
    
    
    
end