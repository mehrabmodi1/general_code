clear all
close all

direc_lists_mat = [...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_alpha1_lowUS.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_alpha1_lowUS_set2.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_alpha1_lowUS_set3.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_alpha1_lowUS_set4.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_gamma1_lowUS_set1.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_alpha1_lowUS_backward_ctrl.xls'}...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\Yoshi_THnull_G1pedc.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_gamma1_Chrim_stim_lifetime.xls'};...
                    {'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_gamma1_lowUS_MB085C.xls'}...
                    ];
                
n_direc_lists = size(direc_lists_mat, 1);

save_path_base = 'E:\Data\Analysed_data\Manual_ROIs\';
fuse_ROIs = 1;          %0-no, 1-yes. This flattens dim3 of ROI_mat in case there is a multi-patch neuron that needs to be considered as a single object.

%looping through all directory lists and all datasets once to save mean frames and again to save manually determined slow x-y  offsets for each trial, as well as determine bad z-drift trials
for do_over = 1:2
    for direc_list_n = 1:n_direc_lists
        list_direc = direc_lists_mat{direc_list_n, 1};
        [del, curr_direc_list] = xlsread(list_direc, 1);
        n_dirs = size(curr_direc_list, 1);
        direc_counter = 0;
        
        %parsing direc list path for name of direc list
        namei = findstr(list_direc, '\');
        namei = namei(end) + 1;
        dir_list_name = (list_direc(namei:(end - 4)));

        %loop to go through all experiment datasets listed in list file
        for direc_counter = 1:n_dirs
            %% House-keeping
            direc = curr_direc_list{direc_counter, 1};
            direc = [direc, '\'];
            remove_small_tifs(direc);
            prev_direc = pwd;
            cd([direc]);
            dataset_namei = findstr(direc, '\20');
            dataset_name = direc((dataset_namei + 1):end);
            save_path = [save_path_base, dataset_name, '\' ];
           
            %checking if dataset has already been analysed
            old_dir = pwd;
            if isdir(save_path) == 1
                cd(save_path);

                tif_list = dir('*.tif');
                if isempty(tif_list) == 0
                    disp([dataset_name, 'already analysed. skipping...'])
                    continue
                else
                end
            else
            end
               
            
            disp(['Reading in avg stack for ', direc])
            
            if do_over == 1
                cd(direc)
                tif_list = dir('*.tif');
                if exist([save_path, '\tr_avg_stack.mat']) == 2
                    ave_stack = load([save_path, '\tr_avg_stack.mat']);
                    ave_stack = ave_stack.ave_stack;
                    if size(ave_stack, 3) < length(tif_list)
                        tif_start_n = size(ave_stack, 3) + 1;
                    else
                        
                        continue
                    end
                else
                    tif_start_n = 1;
                end
                for tif_n = tif_start_n:length(tif_list)
                    curr_stack = ScanImageTiffReader([direc, tif_list(tif_n).name]).data();
                    curr_stack = permute(curr_stack,[2 1 3]);
                    curr_stack = double(curr_stack);
                    ave_stack(:, :, tif_n) = std(curr_stack, 0, 3, 'omitnan');
                    disp(['Saving avg stack, tr ', int2str(tif_n), ' done.'])
                end
                
                mkdir(save_path);
                save([save_path, '\tr_avg_stack.mat'], 'ave_stack');
                clear ave_stack
                
            elseif do_over == 2
                
                dataset_stack = load([save_path, '\tr_avg_stack.mat']);
                dataset_stack = dataset_stack.ave_stack;
                
                %looping through an ave frame for each trial in the
                %dataset, for the user to click on a fixed landmark in each frame
                %to correct x-y drift, and also indicate bad z-trials
                if exist([save_path, '\xy_lags.mat']) == 2
                    continue
                else
                end
                
                done_marking = 0;
                while done_marking == 0
                    [lag_mat, bad_trs, done_marking] = manual_xylags_zbad(dataset_stack);
                end
                
                bad_tr_list = 1:1:size(dataset_stack, 3);
                bad_tr_list(bad_trs == 1) = [];
                save([save_path, '\xy_lags.mat'], 'lag_mat');
                save([save_path, '\bad_trial_list.mat'], 'bad_tr_list');        %bad_tr_list is actually the list of good trials.
                clear lag_mat
                clear bad_trs
                
            end

        end

        
    end
end

%loop to go through all directory lists
for direc_list_n = 1:n_direc_lists
    list_direc = direc_lists_mat{direc_list_n, 1};
    [del, curr_direc_list] = xlsread(list_direc, 1);
    n_dirs = size(curr_direc_list, 1);
    direc_counter = 0;
    
    %parsing direc list path for name of direc list
    namei = findstr(list_direc, '\');
    namei = namei(end) + 1;
    dir_list_name = (list_direc(namei:(end - 4)));
    
    %loop to go through all experiment datasets listed in list file
    for direc_counter = 1:n_dirs
        %% House-keeping
        direc = curr_direc_list{direc_counter, 1};
        direc = [direc, '\'];
        remove_small_tifs(direc);
        prev_direc = pwd;
        cd([direc]);
        tif_list = dir('*.tif');
        
        disp(['Extracting traces from ' direc]);
        try
            curr_stack = ScanImageTiffReader([direc, tif_list(1).name]).data();
        catch
            keyboard
        end
        curr_stack = permute(curr_stack,[2 1 3]);
        ref_im = mean(curr_stack, 3, 'omitnan');
        %loading in manually drawn, FIJI ROIs
        prev_direc = pwd;
        cd([direc, 'ROIs\']);
        ROI_list = dir('*.roi');
        cd(prev_direc);
        n_ROIs = size(ROI_list, 1);
        
        ROI_mat = zeros(size(ref_im, 1), size(ref_im, 2), n_ROIs);
        for ROI_n = 1:n_ROIs
            curr_name = ROI_list(ROI_n).name;
            curr_ROI = ReadImageJROI([direc, 'ROIs\', curr_name]);
            ROI_mat(:, :, ROI_n) = poly2mask(curr_ROI.mnCoordinates(:, 1), curr_ROI.mnCoordinates(:, 2), size(ref_im, 1), size(ref_im, 2));
        end
        
        if fuse_ROIs == 1
            ROI_mat = sum(ROI_mat, 3);
        else
        end
        
        %extracting raw traces
        dataset_namei = findstr(direc, '\20');
        dataset_name = direc((dataset_namei + 1):end);
        save_path = [save_path_base, dataset_name, '\' ];
        
        [raw_data_mat] = extract_raw_traces_par(direc, ROI_mat, save_path, 2);
        
        
        %copying over files needed for further analysis to results
        %directory
        prev_direc = pwd;
        cd([direc])
        %copying over stimulus param file from raw data directory
        param_names = dir('params*.*');
        for p_file_n = 1:size(param_names, 1)
            copyfile([direc, param_names(p_file_n).name], [save_path, param_names(p_file_n).name]);            
        end
        
        % copying over PID traces to the results folder
        PID_fnames = dir('PID*.*');
        for PID_trace_n = 1:size(PID_fnames, 1)
            curr_name = PID_fnames(PID_trace_n).name;
            copyfile([direc, curr_name], [save_path, curr_name]);

        end
        %copying over first trial .tif file to results folder
        tif_fnames = dir('*.tif');
        curr_name = tif_fnames(1).name;
        copyfile([direc, curr_name], [save_path, curr_name]);
        
        disp(['Done extracting traces from ' direc]);
        cd(prev_direc)        
    end
    
end