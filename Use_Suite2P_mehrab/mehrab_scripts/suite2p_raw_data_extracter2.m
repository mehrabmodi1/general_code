 clear all
close all

results_direc = 'E:\Data\Analysed_data\Suite2p\Results\';
reg_tif_direc = 'E:\Data\Analysed_data\Suite2p\Reg_Tiff\';
raw_direc_base = 'E:\Data\Raw_Data_Current\Resonant\';

%% going through raw_direc_base and building a list of all raw data direcs that have a stimulus params file to analyse them
[raw_direc_list] = setup_raw_direc_list(raw_direc_base);

%%Setting up loop to repeatedly run Suite2P for each identified dataset directory (ie. all sub-directories of raw_direc_base that contain stim params files.)
n_direcs_analysed = 0;
for raw_direc_n = 1:size(raw_direc_list, 1)
    raw_direc = raw_direc_list{raw_direc_n, 1};

    %% Pre-prep
    %raw_direc = [raw_direc '\'];
    %Getting rid of empty .tiff files in raw data folder(s) left behind by ScanImage
    remove_small_tifs([raw_direc_base, raw_direc]);

    %setting up master_file.m and make_db.m
    ft_factor = setup_Suite2P_files(raw_direc_base, raw_direc, 'E:\Code\general_code_repo\Use_Suite2P_mehrab\mehrab_scripts\Suite2P_starter_files\master_file.m');

    %Checking if this directory has already been analysed with Suite2P
    if isdir([results_direc raw_direc]) == 0
        prev_direc = pwd;
        cd([raw_direc_base, raw_direc]);
        disp(['currently analysing ' raw_direc_base, raw_direc]);
        
        %running Suite2P
        try
            master_file
        catch
            keyboard
        end
        n_direcs_analysed = n_direcs_analysed + 1;

        cd(prev_direc);
    else
        disp([raw_direc_base, raw_direc ' has already been analysed... skipping.'])
    end
    disp(['analysed ' int2str(n_direcs_analysed) ' new datasets.']);
end


%% loop to repeatedly run manual z-drift detection GUI for all new datasets
for raw_direc_n = 1:size(raw_direc_list, 1)    
    raw_direc = raw_direc_list{raw_direc_n, 1};
    %% Reading in registered tiffs and displaying for manual removal of trials with z-drift
    if exist([results_direc, raw_direc, 'bad_trial_list.mat']) ~= 2 &&  exist([results_direc, raw_direc, 'good_trial_list.mat']) ~= 2
        %[bad_tr_list] = find_bad_trials_res([reg_tif_direc, raw_direc 'Plane1\']);  %these are actually the good trials
        [good_tr_list, bad_tr_list] = find_bad_trials_res([raw_direc_base, raw_direc]);  %these are actually the good trials
        save([results_direc, raw_direc, 'good_trial_list.mat'], 'good_tr_list');
    else
        disp('z-drift trials have already been manually identified... skipping.')
    end
end


%% Loop to repeatedly call ROI_prune to manually curate ROIs
for raw_direc_n = 1:size(raw_direc_list, 1)
    raw_direc = raw_direc_list{raw_direc_n, 1};
    disp(raw_direc);

    %% Loading in Suite2P analysis results file
    prev_direc = pwd;
    cd([results_direc, raw_direc])
    
    if exist([results_direc, raw_direc, '\ROIs_pruned.txt']) ~= 2
        new_main       %this is the ROI pruning GUI that comes with Suite2P
        keyboard
        del = [];
        save([results_direc, raw_direc, '\ROIs_pruned.txt'], 'del');
    else
        disp([results_direc, raw_direc, ' has ROIs already pruned. Skipping.']);
    end
    
    %checking if directory structure was extended with a \1\ folder at the
    %end and copying over the results file from ROI_prune to that folder
    copy_files_to_1_direc(results_direc, raw_direc);
    
end

log_m = [];
%% Extracting raw F data, writing to file
for raw_direc_n = 1:size(raw_direc_list, 1)
    raw_direc = raw_direc_list{raw_direc_n, 1};
    disp(raw_direc);
    cd([results_direc, raw_direc])
    
    %checking if data has already been extracted
    if exist([results_direc, raw_direc, 'trace_extraction_complete.mat']) == 2
       
       disp(['traces have alredy been extracted for ' [results_direc, raw_direc] '. Skipping...']);
       continue
       
    else
        
    end
    
    %reading in most recent Suite2P results file to get Suite2P ROIs
    newest_results_file = find_newest_file([results_direc, raw_direc], '_proc');
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

    %% Extracting raw fluorescence traces after doing a slow xy-correction, and copying over files needed for further analysis
    prev_direc = pwd;
    [raw_data_mat] = extract_raw_traces([raw_direc_base, raw_direc], ROI_mat, [results_direc, raw_direc, '\'], 0);
    disp(['done extracting traces from ' reg_tif_direc, raw_direc]);
    
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