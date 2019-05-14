clear all
close all

results_direc = 'E:\Data\Analysed_data\Suite2p\Results\';
results_direc_manual_ROIs = 'E:\Data\Analysed_data\Manual_ROIs\';
reg_tif_direc = 'E:\Data\Analysed_data\Suite2p\Reg_Tiff\';
raw_direc_base = 'E:\Data\Raw_Data_Current\Resonant\';

extract_both_channels = 0;

%% Reading in manually created direc list
direc_list = 'E:\Data\Raw_Data_Current\dataset_lists\KC_Ca_alpha1T_set2.xls';
[del, raw_direc_list] = xlsread(direc_list, 1);
raw_direc_list_copy = raw_direc_list;

%% Setting up loop to repeatedly run Suite2P for each identified dataset directory (ie. all sub-directories of raw_direc_base that contain stim params files.)
n_direcs_analysed = 0;
rem_dir_list = [];
for raw_direc_n = 1:size(raw_direc_list, 1)
   
    % House-keeping
    direc = raw_direc_list{raw_direc_n, 1};
    direc = [direc, '\'];
    remove_small_tifs(direc);
    prev_direc = pwd;
    cd([direc]);
    dataset_namei = findstr(direc, '\20');
    raw_direc = direc((dataset_namei + 1):end);
    raw_direc_list{raw_direc_n, 1} = raw_direc;    
    %% Pre-prep
    %raw_direc = [raw_direc '\'];
    
    %Checking if this directory has already been analysed with Suite2P or ManualROIextracter

    if isdir([results_direc_manual_ROIs, raw_direc]) == 1
        %raw_direc_list(raw_direc_n) = [];
        continue
    else
    end
    
    if isdir([results_direc raw_direc]) == 0
        prev_direc = pwd;
        cd([raw_direc_base, raw_direc]);
        disp(['currently analysing ' raw_direc_base, raw_direc]);
        
        %Getting rid of empty .tiff files in raw data folder(s) left behind by ScanImage
        remove_small_tifs([raw_direc_base, raw_direc]);

        %setting up directory structure for Suite2P
        raw_direc = setup_raw_direc(raw_direc_base, raw_direc);
        raw_direc_list{raw_direc_n, 1} = raw_direc;
        
        %skipping if current directory doesn't meet criteria in setup_raw_direc
        if isempty(raw_direc) == 1
            continue
        else
        end
        
        %setting up master_file.m and make_db.m
        ft_factor = setup_Suite2P_files(raw_direc_base, raw_direc, 'E:\Code\general_code_repo\Use_Suite2P_mehrab\mehrab_scripts\Suite2P_starter_files\master_file.m');
        
        prev_dir = pwd;
        cd([raw_direc_base, raw_direc]);
        %running Suite2P
%         try
            master_file
%         catch
%             keyboard
%         end
        n_direcs_analysed = n_direcs_analysed + 1;

        cd(prev_direc);
    else
        disp([raw_direc_base, raw_direc ' has already been analysed... skipping.'])
        rem_dir_list = [rem_dir_list; raw_direc_n];
    end
    disp(['analysed ' int2str(n_direcs_analysed) ' new datasets.']);
end
%raw_direc_list(rem_dir_list) = [];



%% Loop to repeatedly call ROI_prune to manually curate ROIs identified by Suite2P
for raw_direc_n = 1:size(raw_direc_list, 1)
    raw_direc = raw_direc_list{raw_direc_n, 1};
    raw_direc = raw_direc_with_1(raw_direc_base, raw_direc);
    disp(raw_direc);

    %% Loading in Suite2P analysis results file
    prev_direc = pwd;
    cd([results_direc, raw_direc])
    
    if exist([results_direc, raw_direc, '\ROIs_pruned.txt']) ~= 2
        new_main       %this is the ROI pruning GUI that comes with Suite2P
        keyboard       %don't remove this keyboard - it's needed for normal running.
        del = [];
        save([results_direc, raw_direc, '\ROIs_pruned.txt'], 'del');
    else
        disp([results_direc, raw_direc, ' has ROIs already pruned. Skipping.']);
    end
    
    %checking if directory structure was extended with a \1\ folder at the
    %end and copying over the results file from ROI_prune to that folder
    copy_files_to_1_direc(results_direc, raw_direc);
    
end



%% Doing a slow, manual xy correction and z-drift detection
%looping through all datasets once to save mean frames and again to save manually determined slow x-y  offsets for each trial, as well as determine bad z-drift trials

for do_over = 1:2
    
    %loop to go through all experiment datasets listed in list file
    for direc_counter = 1:size(raw_direc_list, 1)
        if exist([raw_direc_base, raw_direc, 'skip_direc.txt']) == 2
            continue
        else
        end
        
        %% House-keeping
        raw_direc = raw_direc_list{direc_counter, 1};
        raw_direc = [raw_direc_base, raw_direc];
        remove_small_tifs(raw_direc);
        prev_direc = pwd;
        cd([raw_direc]);
        dataset_namei = findstr(raw_direc, '\20');
        dataset_name = raw_direc((dataset_namei + 9):end);
        raw_direc = raw_direc((dataset_namei + 1):end);
        raw_direc = raw_direc_with_1(raw_direc_base, raw_direc);
        raw_direc = [raw_direc, '\'];
        save_path = [results_direc, raw_direc, '\' ];
        
        %checking if dataset has already been analysed
        old_dir = pwd;
        if isdir(save_path) == 1
            cd(save_path);
                        
            tif_list = dir('*.tif');
            
            %using criterion that if tr 1 has been copied over, this
            %dataset has been aalysed already.
            if isempty(tif_list) == 0
                disp([dataset_name, 'already analysed. skipping...'])
                continue
            else
            end
        else
        end


        disp(['Reading in avg stack for ', raw_direc])

        if do_over == 1
            cd([raw_direc_base, raw_direc, '\']);
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
                stack_obj = ScanImageTiffReader([raw_direc_base, raw_direc, '\', tif_list(tif_n).name]);
                curr_stack = stack_obj.data();
                curr_stack = permute(curr_stack,[2 1 3]);
                curr_stack = double(curr_stack);
                
                [frame_time, zoom, n_chans] = SI_tif_info(stack_obj);
                if n_chans == 2
                    red_stack = curr_stack(:, :, [2:2:end]);
                    curr_stack = curr_stack(:, :, [1:2:end]);

                else
                end
                    
                
                
                
                ave_stack(:, :, tif_n) = std(curr_stack, 0, 3, 'omitnan');
                disp(['Saving avg stack, tr ', int2str(tif_n), ' done.'])
            end
            mkdir(save_path);
            save([save_path, '\tr_avg_stack.mat'], 'ave_stack');
            clear ave_stack

        elseif do_over == 2
            
            %reading in most recent Suite2P results file to get Suite2P ROIs
            newest_results_file = find_newest_file([results_direc, raw_direc], '_proc');
            if isempty(newest_results_file) == 1
                continue
            else
            end
            disp(['Loading Suite2P results file ' newest_results_file])
            data_mat = load([results_direc, raw_direc, newest_results_file]);

            try
                data_mat = data_mat.dat;
                disp('Done loading.')
            catch
                del = [];
                save([raw_direc_base, raw_direc, 'skip_direc.txt'], 'del');
                disp(['no data in ' results_direc, raw_direc, dir_contents(max_datenum(2)).name, '... skipping.']);
                continue
            end
            cd(prev_direc);
            
            %Creating ROI matrix that includes only ROIs manually classified as real cells
            [ROI_mat] = setup_Suite2P_ROIs(data_mat);
            clear data_mat

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
            ROI_mat_s = sum(ROI_mat, 3);
            ROI_mat_s(ROI_mat_s > 0) = 1;
            while done_marking == 0
                [lag_mat, bad_trs, done_marking, bk_ROI] = manual_xylags_zbad2(dataset_stack, ROI_mat_s);
            end

            bad_tr_list = 1:1:size(dataset_stack, 3);
            bad_tr_list(bad_trs == 1) = [];
            save([save_path, '\xy_lags.mat'], 'lag_mat');
            save([save_path, '\bad_trial_list.mat'], 'bad_tr_list');        %bad_tr_list is actually the list of good trials.
            save([save_path, '\bk_ROI.mat'], 'bk_ROI'); 
            save([save_path, '\ROI_mat.mat'], 'ROI_mat');
            clear lag_mat
            clear bad_trs
            
        end

    
    end
end


log_m = [];

%PICK UP THREAD HERE
%fix raw_direc path issue here. figure out why program is moving directly
%to fly3.

%% Extracting raw fluorescence traces after doing a slow xy-correction, and copying over files needed for further analysis
for raw_direc_n = 1:size(raw_direc_list, 1)
    raw_direc = raw_direc_list{raw_direc_n, 1};
    raw_direc = [raw_direc_base, raw_direc];
    remove_small_tifs(raw_direc);
    prev_direc = pwd;
    cd([raw_direc]);
    dataset_namei = findstr(raw_direc, '\20');
    dataset_name = raw_direc((dataset_namei + 1):end);
    raw_direc = raw_direc((dataset_namei + 1):end);
    raw_direc = raw_direc_with_1(raw_direc_base, raw_direc);
    raw_direc = [raw_direc, '\'];
    save_path = [results_direc, raw_direc, '\' ];
    direc = [raw_direc_base, raw_direc];
    
    prev_direc = pwd;
    cd([direc]);
    tif_list = dir('*.tif');

    if exist([results_direc, raw_direc, '\trace_extraction_complete.mat']) == 2
       continue

    else
    end

    disp(['Extracting traces from ' direc]);
    
    %loading in previously saved ROI_mat
    ROI_mat = load([results_direc, raw_direc, '\ROI_mat.mat']);
    ROI_mat = ROI_mat.ROI_mat;


    %extracting raw traces
    dataset_namei = findstr(direc, '\20');
    dataset_name = direc((dataset_namei + 1):end);
    save_path = [results_direc, raw_direc, '\' ];
  

    [raw_data_mat] = extract_raw_traces_par(direc, ROI_mat, save_path, 1, extract_both_channels);
    
    %copying over stimulus param files
    cd([raw_direc_base, raw_direc]);
    param_direc = dir('params*.*');
    for param_file_n = 1:size(param_direc, 1)
        copyfile([raw_direc_base, raw_direc, '\', param_direc(param_file_n).name], [results_direc, raw_direc, '\', param_direc(param_file_n).name]);
    end
    log_m = [log_m; {['finished extracting data for ' raw_direc]}];
    save([results_direc, '\log_file.mat'], 'log_m');
    cd(prev_direc)
    
    %copying over PID traces to the results folder
    cd([raw_direc_base, raw_direc])
    PID_fnames = dir('PID*.*');
    for PID_trace_n = 1:size(PID_fnames, 1)
        curr_name = PID_fnames(PID_trace_n).name;
        copyfile([raw_direc_base, raw_direc, curr_name], [results_direc, raw_direc, curr_name]);
    end
    
    %copying over first trial .tif file to results folder
    cd([raw_direc_base, raw_direc])
    tif_fnames = dir('*.tif');
    curr_name = tif_fnames(1).name;
    copyfile([raw_direc_base, raw_direc, curr_name], [results_direc, raw_direc, curr_name]);

    
    
end