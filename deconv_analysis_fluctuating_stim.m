clear all
close all

direc_lists_mat =  [{'D:\Data\Janelia\resonant\dataset_list_fluc_stim_somatic_20171220.xls'}...
                    
                   ]; 
                
            
n_direc_lists = size(direc_lists_mat, 1);
                
                
color_vec = load('C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Code\std_color_vec.txt');
a = colormap('bone');
greymap = flipud(a);
colormap(greymap)
suppress_plots = 0;       %1 - doesn't plot stuff, 0 - plots stuff

[del, odor_names] = xlsread('D:\Data\CSHL\odor_names_20161108.xls', 1);

%loop to go through all directory lists
for direc_list_n = 1:n_direc_lists

    list_direc = direc_lists_mat{direc_list_n, 1};
    [del, curr_direc_list] = xlsread(list_direc, 1);
    n_dirs = size(curr_direc_list, 1);
    direc_counter = 0;
    
    %parsing direc list path for name of direc list
    namei = findstr(list_direc, 'dataset_list');
    dir_list_name = (list_direc(namei:(end-4)));
        
    %loop to go through all experiment datasets listed in list file
    for direc_counter = 1:n_dirs
        direc = curr_direc_list{direc_counter, 1};
        
        direc = [direc, '\'];
        
        %loading in and parsing params file to get stimulus parameter
        %details
        tif_times = load([direc, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
        tif_times = tif_times.time_stamps;
        [stim_mat, stim_mat_simple, column_heads] = load_params_trains(direc, tif_times);
        odor_list = unique(stim_mat_simple(:, 2) );
        n_odors = length(odor_list);
        odor_dur_list = unique(stim_mat_simple(:, 3) );
        n_od_durs = length(odor_dur_list);
                        
        %loading Suite2P results file
        Suite2P_results = load_suite2P_results(direc);
        frame_rate = Suite2P_results.ops.imageRate;
        frame_time = 1./frame_rate.*1000;     %in ms
        stim_time = stim_mat_simple(1, 7);
        stim_fr = round((stim_time.*1000)./frame_time);
        post_od_scan_dur = stim_mat_simple(1, 10);
        %building a matrix where each row has [stim_fr, stim_end_fr,
        %post_od_scan_end_fr] with one such row for each odor duration.
        stim_frs = [];
        for dur_n = 1:n_od_durs
            curr_dur = odor_dur_list(dur_n);
            stim_frs = [stim_frs; [stim_fr, (round( ((stim_time + curr_dur).*1000)./frame_time)),...
                round( ((stim_time + curr_dur + post_od_scan_dur).*1000)./frame_time)]];
        end
        keyboard        
        %loading extracted raw fluorescence data matrices written by
        %raw_dff_extractor
        raw_data_mat = load([direc 'extracted_raw_data_mat.mat']);
        raw_data_mat = raw_data_mat.raw_data_mat;           %raw F traces extracted from ROIs
        raw_data_mat = raw_data_mat(:, :, stim_mat_simple(:, 1));       %making sure only time-stamp matched trials are used for further analysis
        
        %calculating dF/F traces from raw data
        filt_time = 200;            %in ms, the time window for boxcar filter for generating filtered traces
        [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, direc);
        
        
        keyboard
        
        
        
    end
end
       