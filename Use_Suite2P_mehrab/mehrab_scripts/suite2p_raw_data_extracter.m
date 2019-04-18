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
        master_file
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
    if exist([results_direc, raw_direc, 'bad_trial_list.mat']) ~= 2
        [bad_tr_list] = find_bad_trials_res([reg_tif_direc, raw_direc 'Plane1\']);  %these are actually the good trials
        save([results_direc, raw_direc, 'bad_trial_list.mat'], 'bad_tr_list');
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
%% Extracting raw F data, writing to file, downsampling registered .tiffs in time
for raw_direc_n = 1:size(raw_direc_list, 1)
    raw_direc = raw_direc_list{raw_direc_n, 1};
    disp(raw_direc);
    cd([results_direc, raw_direc])
    
    %checking if data has already been extracted
    if exist([results_direc, raw_direc, 'extracted_raw_data_mat.mat']) == 2
       dir_cont = dir([results_direc, raw_direc, 'extracted_raw_data_mat.mat']);
       if dir_cont.bytes > 1e6
           disp(['traces have alredy been extracted for ' [results_direc, raw_direc] '. Skipping...']);
           continue
       elseif exist([results_direc, raw_direc, 'extraction_completed.txt']) == 2
           disp(['traces have alredy been extracted for ' [results_direc, raw_direc] '. Skipping...']);
           continue
       end
      
    else
    end
    
    dir_contents = dir;
    dir_contents(1:2) = [];
    n_files = size(dir_contents, 1);
    %finding most recent file
    max_datenum = [];
    for file_n = 1:n_files
        fname = dir_contents(file_n).name;
        if isempty(findstr(fname, '_proc')) == 1
            continue
        else
        end
        curr_datenum = dir_contents(file_n).datenum;
        if isempty(max_datenum) == 1
            max_datenum = [curr_datenum, file_n];
        elseif isempty(max_datenum) == 0
            if max_datenum(1) < curr_datenum == 1
                max_datenum(1) = curr_datenum;
                max_datenum(2) = file_n;
            else
            end
        end

    end
    
    try
        disp(['Loading Suite2P results file ' dir_contents(max_datenum(2)).name])
    catch
        keyboard
    end
    
    data_mat = load([results_direc, raw_direc, dir_contents(max_datenum(2)).name]);
    try
        data_mat = data_mat.dat;
    catch
        del = [];
        save([raw_direc_base, raw_direc, 'skip_direc.txt'], 'del');
        disp(['no data in ' results_direc, raw_direc, dir_contents(max_datenum(2)).name, '... skipping.']);
        continue
    end
    disp('Done loading.')

    cd(prev_direc);

    %% Creating ROI matrix that includes only ROIs manually classified as real cells
    [ROI_mat] = setup_Suite2P_ROIs(data_mat);
    clear data_mat

    %% Extracting raw fluorescence traces from registered tiff images
    prev_direc = pwd;
    
    disp('Extracting traces from raw_direc rather than reg_tif_direc.')
    extraction_type = 1;        %0 - extract traces from reg .tifs; 1 - extract traces from raw .tifs    
    
    if extraction_type == 0
        cd([reg_tif_direc, raw_direc, '\Plane1\']);
        t_num_template = '_plane';
        t_num_gap = 8;
    elseif extraction_type == 1
        cd([raw_direc_base, raw_direc]);
        t_num_template = '_00';
        t_num_gap = 2;
    else
    end
        
    dir_contents = dir('*.tif');
    n_trials = size(dir_contents, 1);
    t_num_list = zeros(n_trials, 1) + nan;
    
    %preparing to extract raw F data
    n_frames_list = zeros(n_trials, 1) + nan;
    
    for trial_n = 1:n_trials
        curr_name = dir_contents(trial_n).name;
        try
            t_numi = findstr(curr_name, t_num_template);
            t_doti = findstr(curr_name, '.tif');
        catch
            error(['Error parsing registered tiff filename to find trial number. Search string ', t_num_template, ' not found in current filename.'])
        end
        t_num = str2num(curr_name( (t_numi + t_num_gap):(t_doti - 1) ));        %ScanImage acquisition number for current trial filename
        t_num_list(trial_n, 1) = t_num;

        try
            stack_obj = ScanImageTiffReader([reg_tif_direc, raw_direc, curr_name(5:end)]);
            desc = stack_obj.descriptions;
            n_frames_list(trial_n, 1) = size(desc, 1);
        catch
            n_frames_list(trial_n, 1) = 0;
            
        end
        
        clear stack_obj;
        clear desc;
        disp(['building list of frame_ns; trial ', int2str(trial_n), ' of ', int2str(n_trials), ' done.']);
    end

    n_frames_max = max(n_frames_list);                    %making sure there's room in the storage matrix for all he frames from the longest trials

    t_num_list = sort(t_num_list);                    %sorted list of ScanImage acq numbers - used to make sure trials are read in in the correct order
    n_ROIs = size(ROI_mat, 3);
    raw_data_mat = zeros(n_frames_max, n_ROIs, n_trials) + nan;
    
    
    
    for trial_n = 1:n_trials
        curr_t_num = t_num_list(trial_n);       

        %identifying the trial name that should be read in next
        warning = 1;
        for trial_n_sub = 1:n_trials
            curr_name = dir_contents(trial_n_sub).name;
            t_numi = findstr(curr_name, t_num_template);
            t_doti = findstr(curr_name, '.tif');
            t_num = str2num(curr_name( (t_numi + t_num_gap):(t_doti - 1) ));

            if t_num == curr_t_num
                warning = 0;
                break
            else
            end

        end
        
        if warning == 1
            error('Could not find correct ScanImage file number to read in.');
        else
        end

        if extraction_type == 0
            stack = load_all_frames([reg_tif_direc, raw_direc, 'Plane1\' curr_name]);
        elseif extraction_type == 1
            stack = load_all_frames([raw_direc_base, raw_direc, curr_name]);
        else
        end
        n_frames = size(stack, 3);
        log_m = [log_m; {['loaded stack ' curr_name]}];
        save([results_direc, '\log_file.mat'], 'log_m');
    %     m_data=ScanImageTiffReader([reg_tif_direc, raw_direc, 'Plane1\' curr_name]);
    %     m_data = m_data.metadata;

        %reading in raw fluorescence data for each ROI into a single matrix (dim1 - frame_n, dim2 - ROI_n)
        raw_vec = zeros(1, n_ROIs);
        for frame_n = 1:n_frames
            curr_frame = stack(:, :, frame_n);
            curr_frame = im2double(curr_frame);

            for ROI_n = 1:n_ROIs
                curr_ROI = ROI_mat(:, :, ROI_n);
                raw_vec(1, ROI_n) = sum(sum(curr_frame.*curr_ROI));

            end
            raw_data_mat(frame_n, :, trial_n) = raw_vec;
        end
        clear stack
        disp(['traces extracted from trial ', int2str(trial_n)])
        log_m = [log_m; {['extracted traces from ' curr_name]}];
        save([results_direc, '\log_file.mat'], 'log_m');
    end
    
    %Saving raw_data_mat to disk and copying the stimulus params file over to
    %the results directory - all set for further analysis after this.
    save([results_direc, raw_direc, '\extracted_raw_data_mat.mat'], 'raw_data_mat');
    disp(['done extracting traces from ' reg_tif_direc, raw_direc]);
    cd([raw_direc_base, raw_direc]);
    param_direc = dir('params*.*');
    for param_file_n = 1:size(param_direc, 1)
        copyfile([raw_direc_base, raw_direc, '\', param_direc(param_file_n).name], [results_direc, raw_direc, '\', param_direc(param_file_n).name]);
    end
    log_m = [log_m; {['finished extracting data for ' raw_direc]}];
    save([results_direc, '\log_file.mat'], 'log_m');
    cd(prev_direc)
    
    %saving timestamps for each tif file in the current dataset
    time_stamps = read_tiff_timestamps([raw_direc_base, raw_direc]);
    save([results_direc, raw_direc, '\tif_time_stamps.mat'], 'time_stamps');
    
    
    %% downsampling registered tiffs in time and over-writing old ones. deleting raw, registered frames.
    %tiff_downsampler([reg_tif_direc, raw_direc], round(10./ft_factor));
    
    %% copying over PID traces to the results folder
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

    del = [];
    save([results_direc, raw_direc, '\extraction_completed.txt'], 'del');       %saving a dummy file to indicate that all the raw data was extracted successfully.
    
end