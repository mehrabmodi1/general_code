clear all
close all

direc_lists_mat =  [{'C:\Data\Data\Analysed_data\dataset_list_fluc_stim_somas_20171226.xls'}...
                    
                   ]; 
                
            
n_direc_lists = size(direc_lists_mat, 1);
                
global color_vec;                
color_vec = load('C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Code\std_color_vec.txt');
a = colormap('bone');
global greymap
greymap = flipud(a);
colormap(greymap)
suppress_plots = 0;       %1 - doesn't plot stuff, 0 - plots stuff

kernel_width = 20;         %in s, the duration of dF/F trace to use and extract as the Ca-response kernel.

[del, odor_names] = xlsread('C:\Data\Data\Analysed_data\odor_names_20161108.xls', 1);

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
        
        %identifying significantly responsive cells
        [resp_areas, sig_trace_mat, sig_trace_mat_old, sig_cell_mat] = cal_sig_responses_res(dff_data_mat, stim_mat, stim_mat_simple, direc, frame_time);
        
        %Running data quality control checks
        sig_cell_mat_old = sig_cell_mat;
        [sig_cell_mat, all_bad_trs] = cell_data_quality_control(dff_data_mat_f, stim_mat, stim_mat_simple, sig_cell_mat, 0);
        
        
        %% Running fitting algorithm
        for odor_n = 1:n_odors
            odor_ni = odor_list(odor_n);
            curr_sig_cells = find(sig_cell_mat(:, odor_ni) == 1);
            for train_n = 1:n_trains
                curr_trs = find(stim_mat_simple(:, 2) == odor_ni & stim_mat_simple(:, 11) == train_n);
                ave_dff_resp_mat = mean(dff_data_mat_f(:, curr_sig_cells, curr_trs), 3, 'omitnan');
                %normalising each dff response trace
                max_vec = max(ave_dff_resp_mat, [], 1);
                max_mat = repmat(max_vec, size(ave_dff_resp_mat, 1), 1);
                ave_dff_resp_mat = ave_dff_resp_mat./max_mat;
                
                PID_traces = get_PID_traces(direc, curr_trs, frame_time);        %getting mean PID trace
                PID_trace_mean = mean(PID_traces, 1)';
                
                %matching sizes of Ca-resps and mean PID trace
                length_diff = size(PID_trace_mean, 1) - size(ave_dff_resp_mat, 1);
                if sign(length_diff) == 1       %ie PID trace is longer
                    PID_trace_mean = PID_trace_mean(1:(end - length_diff), :);
                elseif sign(length_diff) == -1  %ie Ca traces are longer
                    ave_dff_resp_mat = ave_dff_resp_mat(1:(end - abs(length_diff)), :);
                else
                end
                
                %Identifying shortest stim pulse in the train and using responses 
                %to it as the initial kernel
                curr_train = stim_mat(curr_trs(1)).rand_trains;
                [del, short_pulsei] = min(curr_train(:, 2));
                pulse_frames = compute_pulse_frames_train(curr_train, frame_time, stim_mat(curr_trs(1)).stim_latency);
                pulse_frames = pulse_frames(short_pulsei, :);
                kernel_width_f = kernel_width./(frame_time./1000);      %kernel width in frames
                kernel_pad_l = kernel_width_f - (pulse_frames(2) - pulse_frames(1) );
                kernel_traces = ave_dff_resp_mat(pulse_frames(1):(pulse_frames(2) + kernel_pad_l), :);
                
                for cell_n = 1:length(curr_sig_cells)
                    starter_kernel = kernel_traces(:, cell_n);
                    ave_resp_trace = ave_dff_resp_mat(:, cell_n);
                    fitted_kernel = GetKernel_mehrab(starter_kernel, ave_resp_trace, PID_trace_mean);
                    
                    %calculating predicted ave_resp_trace
                    dff_trace_predic = conv(fitted_kernel, PID_trace_mean);
                    dff_trace_predic_zero = conv(starter_kernel, PID_trace_mean);
                    figure(1)
                    subplot(2, 1, 1)
                    plot(ave_resp_trace, 'b')
                    hold on
                    plot(dff_trace_predic, 'r')
                    %plot(dff_trace_predic_zero, 'g')       %predicted train response base on starter kernel
                    
                    subplot(2, 1, 2)
                    plot(starter_kernel, 'b')
                    hold on
                    plot(fitted_kernel, 'r')
                    del = input('press enter');
                    close all
                    
                end
            end
        end
        keyboard
        
        
        
    end
end
       