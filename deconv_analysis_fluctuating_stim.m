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
suppress_plots = 0;       %0 - doesn't plot quality control stuff, 1 - plots stuff

kernel_width = 5;         %in s, the duration of dF/F trace to use and extract as the Ca-response kernel.
n_trs_to_analyse = 15;

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
    
    MSE_mat = [];
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
        [sig_cell_mat, all_bad_trs] = cell_data_quality_control(dff_data_mat_f, stim_mat, stim_mat_simple, sig_cell_mat, suppress_plots);
        disp([num2str((length(all_bad_trs)./size(dff_data_mat, 3)).*100) ' percent of trials were auto-identified as bad and removed.']);
        dff_data_mat(:, :, all_bad_trs) = nan;
        
        %% Running fitting algorithm
        kernel_mat = [];
        for odor_n = 1:n_odors
            odor_ni = odor_list(odor_n);
            curr_sig_cells = find(sig_cell_mat(:, odor_ni) == 1);
            for train_n = 1:n_trains
                other_train_n = 1:n_trains;
                del = find(other_train_n == train_n);
                other_train_n(del) = [];
                other_train_n = other_train_n(randperm(length(other_train_n)));
                other_train_n = other_train_n(1);        %a randomly selected train number other than the one currently being used for the fit
                curr_trs = find(stim_mat_simple(:, 2) == odor_ni & stim_mat_simple(:, 11) == train_n);
                curr_trs = curr_trs(1:n_trs_to_analyse);
                curr_trs2 = find(stim_mat_simple(:, 2) == odor_ni & stim_mat_simple(:, 11) == other_train_n);
                curr_trs2 = curr_trs2(1:n_trs_to_analyse);                
                ave_dff_resp_mat = mean(dff_data_mat_f(:, curr_sig_cells, curr_trs), 3, 'omitnan');
                ave_dff_resp_mat2 = mean(dff_data_mat_f(:, curr_sig_cells, curr_trs2), 3, 'omitnan');
                %normalising each dff response trace
                max_vec = max(ave_dff_resp_mat, [], 1);
                max_mat = repmat(max_vec, size(ave_dff_resp_mat, 1), 1);
                ave_dff_resp_mat = ave_dff_resp_mat./max_mat;
                max_vec2 = max(ave_dff_resp_mat2, [], 1);
                max_mat2 = repmat(max_vec2, size(ave_dff_resp_mat2, 1), 1);
                ave_dff_resp_mat2 = ave_dff_resp_mat2./max_mat2;
                
                PID_traces = get_PID_traces(direc, curr_trs, frame_time);        %getting PID traces
                PID_traces2 = get_PID_traces(direc, curr_trs2, frame_time);      %getting PID traces for a different, randomly selected train 
                PID_trace_mean = mean(PID_traces, 1)';
                PID_trace_mean2 = mean(PID_traces2, 1)';
                
                
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
                pulse_frames_s = pulse_frames(short_pulsei, :);
                kernel_width_f = kernel_width./(frame_time./1000);      %kernel width in frames
                kernel_pad_l = kernel_width_f - (pulse_frames_s(2) - pulse_frames_s(1) );
                kernel_traces = ave_dff_resp_mat(pulse_frames_s(1):(pulse_frames_s(2) + kernel_pad_l), :);
                
                for cell_n = 1:length(curr_sig_cells)
                    starter_kernel = kernel_traces(:, cell_n);
                    ave_resp_trace = ave_dff_resp_mat(:, cell_n);
                    fitted_kernel = GetKernel_mehrab(starter_kernel, ave_resp_trace, PID_trace_mean);
                                        
                    %calculating predicted ave_resp_trace
                    dff_trace_predic = conv(fitted_kernel, PID_trace_mean);
                    dff_trace_predic = dff_trace_predic(1:size(ave_resp_trace, 1));
                    dff_trace_predic_zero = conv(starter_kernel, PID_trace_mean);
                    
                    %calculating predicted response to other stimulus train with the same fitted kernel
                    ave_resp_trace2 = ave_dff_resp_mat2(:, cell_n);
                    dff_trace_predic2 = conv(fitted_kernel, PID_trace_mean2);
                    dff_trace_predic2 = dff_trace_predic2(1:size(ave_resp_trace2, 1));
                        
                    plotting = 0;    
                    if plotting == 1
                        %plotting fit results w.r.t fitted pulse train
                        figure(1)
                        subplot(3, 1, 1)
                        plot(ave_resp_trace, 'b')
                        hold on
                        plot(dff_trace_predic, 'r')
                        add_stim_bar(1, pulse_frames, color_vec(odor_ni, :))
                        %plot(dff_trace_predic_zero, 'g')       %predicted train response base on starter kernel

                        subplot(3, 1, 2)
                        %plot(starter_kernel, 'b')
                        hold on
                        plot(fitted_kernel, 'r')

                        subplot(3, 1, 3)
                        %plotting fit results w.r.t different pulse train
                        plot(ave_resp_trace2, 'b')
                        hold on
                        plot(dff_trace_predic2, 'r')
                        title('other train convolved with kernel')
                        %plot(dff_trace_predic_zero, 'g')       %predicted train response base on starter kernel
                    
                        MSE1 = mean((ave_resp_trace - dff_trace_predic).^2);
                       
                        if MSE1 < 0.1
                            MSE1
                            keyboard
                            %del = input('press enter for next cell.');
                        else
                        end
                    else
                    end
                    
                    MSE1 = mean((ave_resp_trace - dff_trace_predic).^2);
                    MSE2 = mean((ave_resp_trace2 - dff_trace_predic2).^2);
                    kernel_mat = [kernel_mat,  fitted_kernel];
                    
                    %quantifying cell's own trial to trial variability to ground-truth goodness of fit
                    resp_mat = squeeze(dff_data_mat(:, curr_sig_cells(cell_n), :));
                    ave_trace = mean(resp_mat, 2, 'omitnan');
                    resp_errors = compute_mean_sq_errors(resp_mat, ave_trace);      %vector of MSEs for each trial, with mean response vector
                    
                    %computing signal to noise in mean response trace used to fit kernel
                    base_std = std(ave_resp_trace(1:(stim_mat(curr_trs(1)).stim_latency./(frame_time./1000)) - 1), 1);
                    pk_sig = max(ave_resp_trace);
                    sig_noise = pk_sig/base_std;
                    
                    MSE_mat = [MSE_mat; MSE1, MSE2, mean(resp_errors, 'omitnan'), sig_noise];
                    close all
                    
                end
                
                
                
                %keyboard
                
            end
        end
        %keyboard
        
        
        
    end
    %Comparing goodness of fit for fitted trace with that for independent trace
    figure(3)
    plot(MSE_mat(:, 3), MSE_mat(:, 1), '.')
    xlabel('mean squared error of single tr traces with mean trace')
    ylabel('mean squared error with fitted response trace')
    make_axes_equal(3, 1)
    fig_wrapup(3)

    figure(4)
    plot(MSE_mat(:, 1), MSE_mat(:, 2), '.')
    xlabel('mean squared error with fitted response trace')
    ylabel('mean squared error with independent response trace')
    make_axes_equal(4, 1)
    fig_wrapup(4)
    
    figure(5)
    plot(MSE_mat(:, 4), MSE_mat(:, 1), '.')
    xlabel('dF/F trace signal to noise')
    ylabel('mean squared error with fitted response trace')
    fig_wrapup(5)
    [cor2, p] = corrcoef([MSE_mat(:, 4), MSE_mat(:, 1)])
    keyboard
end
       