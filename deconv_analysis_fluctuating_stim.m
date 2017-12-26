clear all
close all

direc_lists_mat =  [{'C:\Data\Data\Analysed_data\dataset_list_fluc_stim_somas_20171226.xls'}...
                    
                   ]; 
                
            
n_direc_lists = size(direc_lists_mat, 1);
                
                
color_vec = load('C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Code\std_color_vec.txt');
a = colormap('bone');
greymap = flipud(a);
colormap(greymap)
suppress_plots = 0;       %1 - doesn't plot stuff, 0 - plots stuff

%[del, odor_names] = xlsread('D:\Data\CSHL\odor_names_20161108.xls', 1);

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
        n_trains = max(stim_mat_simple(:, 11));
        
        
        %reading ScanImage meta data from a raw .tif
        old_direc = pwd;
        cd(direc);
        tif_names = dir('*.tif');
        cd(old_direc);
        clear old_direc
        tif_name = tif_names(1).name;
        stack_obj = ScanImageTiffReader([direc, tif_name]);
        %metadata = stack_obj.descriptions;
        metadata = stack_obj.metadata;
        fr_ratei = strfind(metadata, 'scanFrameRate');
        frame_rate = metadata((fr_ratei + 16):(fr_ratei + 20));
        frame_rate = str2num(frame_rate);                       %base frame rate in Hz
        avg_factori = strfind(metadata, 'logAverageFactor');
        avg_factor = metadata((avg_factori + 19):(avg_factori + 20));
        avg_factor = str2num(avg_factor);
        frame_rate = frame_rate./avg_factor;
        frame_time = 1./frame_rate.*1000;     %in ms
        stim_time = stim_mat_simple(1, 7);
        stim_fr = round((stim_time.*1000)./frame_time);
        post_od_scan_dur = stim_mat_simple(1, 10);
        
        %loading extracted raw fluorescence data matrices written by
        %raw_dff_extractor
        raw_data_mat = load([direc 'extracted_raw_data_mat.mat']);
        raw_data_mat = raw_data_mat.raw_data_mat;           %raw F traces extracted from ROIs
        raw_data_mat_orig = raw_data_mat;
        raw_data_mat = raw_data_mat(:, :, stim_mat_simple(:, 1));       %making sure only time-stamp matched trials are used for further analysis
        
        %calculating dF/F traces from raw data
        filt_time = 200;            %in ms, the time window for boxcar filter for generating filtered traces
        [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, direc);
        
        %identifying significantly responsive cells
        [resp_areas, sig_trace_mat, sig_trace_mat_old, sig_cell_mat] = cal_sig_responses_res(dff_data_mat, stim_mat, stim_mat_simple, direc, frame_time);
        
        %plotting trial-averaged response traces for each cell
        for odor_n = 1:n_odors
            odor_ni = odor_list(odor_n);
            for train_n = 1:n_trains
                curr_od_trs = find(stim_mat_simple(:, 2) == odor_ni);
                curr_trn_trs = find(stim_mat_simple(:, 11) == train_n);
                curr_trs = intersect(curr_od_trs, curr_trn_trs);           %current trs
                ave_mat = mean(dff_data_mat(:, :, curr_trs), 3, 'omitnan');
                curr_sig_cells = sig_cell_mat(:, odor_ni);
                curr_sig_cells1 = find(curr_sig_cells == 1);
                curr_sig_cells0 = find(curr_sig_cells == 0);
                figure(1)
                imagesc(ave_mat(:, curr_sig_cells1)', [0, 1.2])
                figure(2)
                imagesc(ave_mat(:, curr_sig_cells0)', [0, 1.2])

            end
        end
        keyboard
        
        
        
    end
end
       