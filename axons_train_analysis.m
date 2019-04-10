clear all
close all

direc_lists_mat =  [...
                      %{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_fluc_stim_somas_20171226.xls'}...
                       {'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_fluc_stim_axons_20180117'}, ...
                   ]; 

save_path = 'C:\Data\Data\Analysed_data\Analysis_results\train_clustering\';
n_direc_lists = size(direc_lists_mat, 1);
                
global color_vec;                
color_vec = load('C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Code\std_color_vec.txt');
a = colormap('bone');
global greymap
greymap = flipud(a);
colormap(greymap)
suppress_plots = 0;       %0 - doesn't plot quality control stuff, 1 - plots stuff

kernel_width = 5;         %in s, the duration of dF/F trace to use and extract as the Ca-response kernel.
n_trs_to_analyse = 15;

[del, odor_names] = xlsread('C:\Data\Code\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);




%loop to go through all directory lists
for direc_list_n = 1:n_direc_lists

    list_direc = direc_lists_mat{direc_list_n, 1};
    [del, curr_direc_list] = xlsread(list_direc, 1);
    n_dirs = size(curr_direc_list, 1);
    direc_counter = 0;
    
    %parsing direc list path for name of direc list
    namei = findstr(list_direc, '\');
    namei = namei(end) + 1;
    dir_list_name = (list_direc(namei:(end-4)));
   
    %loop to go through all experiment datasets listed in list file
    for direc_counter = 1:n_dirs
       if direc_counter <= 2
           continue
       else
       end
        
        %% House-keeping
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
        global frame_time;
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
        
        for odor_n = 1:n_odors
            odor_ni = odor_list(odor_n);
            for train_n = 1:n_trains
                curr_trs = find(stim_mat_simple(:, 2) == odor_ni & stim_mat_simple(:, 11) == train_n);
                ave_resp_mat = nanmean(dff_data_mat(:, :, curr_trs), 3);
                
                %normalising ave_resps
                resp_max = max(ave_resp_mat, [], 1);
                ave_resp_mat = ave_resp_mat./repmat(resp_max, size(ave_resp_mat, 1), 1);
                
                figure(1)
                imagesc(ave_resp_mat', [0, 1])
                colormap(greymap)
                curr_train = stim_mat(curr_trs(1)).rand_trains;
                stim_latency = stim_mat(curr_trs(1)).stim_latency;
                pulse_frs = compute_pulse_frames_train(curr_train, frame_time, stim_latency);
                set_xlabels_time(1, (frame_time./1000), 1)
                fig_wrapup(1)
                add_stim_bar(1, pulse_frs, color_vec(odor_n, :));
                
                figure(2)
                plot_vs = 0:(frame_time./1000):(size(dff_data_mat, 1).*(frame_time./1000));
                plot_vs = plot_vs(1:(end - 1));
                cascade_plot(2, squeeze(dff_data_mat(:, 2, curr_trs)), plot_vs, 2, 1, 1, 1, [0.5, 0.5, 0.5])
                add_stim_bar(2, pulse_frs, color_vec(odor_n, :));
                xlabel('time (s)')
                keyboard
                close figure 1
                close figure 2
                
            end
        
        end
        keyboard
    end
end