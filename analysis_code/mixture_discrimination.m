clear all
close all

direc_lists_mat =  [...
                      {'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_mix_disc_somas_20180226.xls'}...
                        %{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_fluc_stim_axons_20180117'}, ...
                   ]; 

save_path = 'C:\Data\Data\Analysed_data\Analysis_results\train_clustering\';
n_direc_lists = size(direc_lists_mat, 1);
                
% global color_vec;                
% color_vec = load('C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Code\std_color_vec.txt');

a = colormap('bone');
global greymap
greymap = flipud(a);
colormap(greymap)
suppress_plots = 0;       %0 - doesn't plot quality control stuff, 1 - plots stuff
[del, odor_names] = xlsread('C:\Data\Code\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);

make_PC_plot = 0;

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
        %odor_list = [odor_list(1); odor_list(3:6); odor_list(2)];       %re-arranging to put the pure odors at opposite ends of the list
        MCHi = find(odor_list == 4);
        odor_list(MCHi) = [];
        odor_list = [odor_list; 4];
        
        n_odors = length(odor_list);
        odor_dur_list = unique(stim_mat_simple(:, 3) );
        n_od_durs = length(odor_dur_list);
        n_trains = max(stim_mat_simple(:, 11));
        
        color_vec = [(linspace(0, 1, n_odors))', (zeros(n_odors, 1) + 0), (linspace(1, 0, n_odors))'];
        
        
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
        dff_data_mat(:, :, all_bad_trs) = nan;
        %disp([num2str((length(all_bad_trs)./size(dff_data_mat, 3)).*100) ' percent of trials were auto-identified as bad and removed.']);
        sig_cells = find(sum(sum(sig_cell_mat, 3), 2) > 0);         %list of all cells significant for any odor for any duration

        
        if make_PC_plot == 1
            %% Plotting PCA points for each odor at each duration
            
            %computing PC weights using only the 60s duration trials

            %concatenating trial time-series for 60s trials to compute PC
            %weights
            curr_trs = find(stim_mat_simple(:, 3) == 60 & stim_mat_simple(:, 12) == 0);    %finding long dur trials that are not train trials

            big_X = [];
            for curr_tr_n = 1:size(curr_trs, 1)
                curr_tr = curr_trs(curr_tr_n);
                curr_traces = squeeze(dff_data_mat_f(:, sig_cells, curr_tr));
                big_X = [big_X; curr_traces];
            end
            [weights, score] = pca(big_X);

            %computing and plotting first two PC values for each odor, at each odor dur
            for dur_n = 1:3
                curr_dur = odor_dur_list(dur_n);
                stim_frs = compute_pulse_frames_train([0, curr_dur], frame_time, stim_mat_simple(1, 7));
                n_acq_frs = ceil((stim_mat_simple(1, 7) + curr_dur + stim_mat_simple(1, 10) )./(frame_time./1000));
                for odor_n = 1:n_odors
                    odor_ni = odor_list(odor_n);
                    curr_trs = find(stim_mat_simple(:, 2) == odor_ni & stim_mat_simple(:, 3) == curr_dur & stim_mat_simple(:, 12) == 0);
                    resp_mat = dff_data_mat_f(1:n_acq_frs, sig_cells, curr_trs);
                    resp_mat_ave = mean(resp_mat, 3, 'omitnan');
                    %PID_traces = get_PID_traces(direc, curr_trs, frame_time);

                    %calculating PC projections
                    nPCs = 2;
                    PC_resp_mat_ave = zeros(n_acq_frs, nPCs) + nan;
                    PC_resp_mat = zeros(n_acq_frs, nPCs, length(curr_trs)) + nan;
                    for PC_n = 1:nPCs
                        curr_weights = weights(:, PC_n)';
                        weight_mat_ave = repmat(curr_weights, size(resp_mat_ave, 1), 1);
                        weight_mat = repmat(weight_mat_ave, 1, 1, length(curr_trs));
                        PC_resp_mat_ave(:, PC_n) = sum((resp_mat_ave.*weight_mat_ave) , 2);
                        PC_resp_mat(:, PC_n, :) = sum((resp_mat.*weight_mat) , 2);
                    end
                    pk_vals_ave = max(PC_resp_mat_ave(stim_frs(1):stim_frs(2), :));
                    pk_vals = max(PC_resp_mat(stim_frs(1):stim_frs(2), :, :));
                    pk_vals = squeeze(pk_vals);
                    n_real_trs = length(find(isnan(pk_vals(1, :)) == 0));
                    pk_vals_se = std(pk_vals, [], 2, 'omitnan')./sqrt(n_real_trs)./2;

                    curr_color = color_vec(odor_n, :);
                    %plotting pk resps in PC-space
                    figure(1)
                    plot(pk_vals_ave(1), pk_vals_ave(2), 'O', 'MarkerFaceColor', curr_color, 'MarkerEdgeColor', 'w', 'MarkerSize', dur_n.*4)
                    hold on
                    e = errorbar(pk_vals_ave(1), pk_vals_ave(2), pk_vals_se(2), '.');
                    e.Color = curr_color;
                    e.CapSize = 0;
                    h = herrorbar(pk_vals_ave(1), pk_vals_ave(2), pk_vals_se(1), pk_vals_se(1), '.');
                    set(h, 'Color', curr_color);
                    xlabel('response peak on PC1') 
                    ylabel('response peak on PC2')  

                end

            end
            hold off
        else
        end
keyboard
        
        %% Template-matching decoder
        %loop to go through trial and compare it with it's duration-matched templates
        for trial_n = 1:size(dff_data_mat, 3)            
            curr_dur = stim_mat_simple(trial_n, 3);
            curr_od_ni = stim_mat_simple(trial_n, 2);
            stim_frs = compute_pulse_frames_train([0, curr_dur], frame_time, stim_mat_simple(1, 7));
            curr_trace_full = dff_data_mat_f((stim_frs(1) - round(5000./frame_time)):(stim_frs(2) + round(5000./frame_time)), sig_cells, trial_n);
            
            for odor_n_tem = 1:n_odors
                odor_ni_tem = odor_list(odor_n_tem);
                curr_trs_tem = find(stim_mat_simple(:, 2) == odor_ni_tem & stim_mat_simple(:, 3) == curr_dur & stim_mat_simple(:, 12) == 0);
                %making sure currently decoded trial is not included in template calculation
                if odor_ni_tem
                    curr_tr_ni = find(curr_trs_tem == trial_n);
                    curr_trs_tem(curr_tr_ni) = [];
                else
                end
                template_full = mean(dff_data_mat_f((stim_frs(1) - round(5000./frame_time)):(stim_frs(2) + round(5000./frame_time)), sig_cells, curr_trs_tem), 3, 'omitnan');
                
                %computing difference between template and current trial for various durations of traces
                step_length = round(1000./frame_time);          %computing n_frames in a t-step of 1s
                n_t_steps = floor( size(template_full, 1)./step_length);     %number of steps in 
                for step_n = 1:n_t_steps
                    curr_trace = curr_trace_full(1:floor(step_length.*step_n), :);
                    curr_template = template_full(1:floor(step_length.*step_n), :);
                    
                    %computing distance between two matrices using distance formula
                    distance = sqrt(sum(sum(sum( (curr_trace - curr_template).^2 ))));
                    curr_tr_dist = sqrt(sum(sum(sum( (curr_trace).^2 ))));
                    curr_temp_dist = sqrt(sum(sum(sum( (curr_template).^2 ))));
                    
                    %PICK UP THREAD HERE
                    %figure out how meaningful this distance is (includes
                    %all non-responders). measure decoder performance.
                    %consider thresholding out non-responsive parts of
                    %traces (?)
                    
                end
            
            end
            
        end
        
        
        
       
        
        
        
        
        keyboard
    end
    
    
end
