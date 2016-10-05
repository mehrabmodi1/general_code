clear all
close all

list_direc = ['D:\Data\CSHL\dataset_list_PN_GH146_20161002.txt'];
fid = fopen(list_direc);

direc_counter = 0;

%loop to go through all experiment datasets listed in list file
while 1
    direc_counter = direc_counter + 1;
    direc = fgetl(fid);
       
    if ischar(direc) ~= 1
        break
    else
    end
    
    direc = [direc '\'];
    
    cd(direc)


    %importing raw .tiff files and saving 2P objects as .mat files; skips if
    %this has already been done before
    importSItree_m();
  

    %loading in 2P objects
    d=dir('*.mat');
    d(find(strcmp({d.name},'.')))=[];
    d(find(strcmp({d.name},'..')))=[];
    d(find(strcmp({d.name},'expt.mat')))=[];
    d(find(strcmp({d.name},'expt_raw_traces.mat')))=[];
    
    if exist([direc 'expt.mat']) ~= 2
        %loop to prep all datasets for data extraction 
        for set_n = 1:length(d)
            curr_set = d(set_n).name;

            %loading 2P object
            dataset = load(curr_set);
            dataset = dataset.data;
            n_trials = size(dataset, 2);

            %loading params.mat file to add stimulus parameters to 2P object
            sub_dir = [direc curr_set(1:(length(curr_set) - 4) ) '\'];
            cd(sub_dir)
            param_name = dir('params*.*');
            params = load(param_name.name);
            params = params.params;


            %adding params to 2P object
            dataset = addStimParamsSI_m(dataset, params);

             try 
                var = dataset(1).process.regSeries;
             catch
                %motion correction of each trial to the frame average of the first 20
                %frames of the trial
                alignStack(dataset, [1,20], 'verbose', 0)

                %updating info in 2P object
                updateMuStack(dataset)

             end
             
             %saving new, combined mat file and continuing analysis with it
             data = dataset;
             save([direc curr_set], 'data');
             clear data
             cd(direc)
            
        end


        %building a single mat file that references all the trials of this
        %experiment
        dataset = combine_matfiles();
        curr_set = 'expt.mat';
        
        
        %saving new, combined mat file and continuing analysis with it
        data = dataset;
        save([direc curr_set], 'data');
        clear data

        %aligning the aligned trials to each other
        alignRepeats(dataset)

        %updating info in 2P object
        updateMuStack(dataset)


        %saving changes to the 2P object
        data = dataset;
        save([direc curr_set], 'data');
        clear data
        
        %adding empty background ROIs to all trials
        ROI.level = 0.023;
        ROI.roi = zeros(255, 256) + 1;
        del_roi = abs(ROI.roi - 1);
        ROI.roi(1:10, 1:10) = 0;                  %tiny patch in top left corner to use as background for calculations
        ROI.notes = 'empty background ROI. DO NOT DELETE ME';
        for tr_n = 1:n_trials
            del = dataset(tr_n).imageStack;
            ROI.backgroundLevel = mean(mean(del(:, :, 10).*del_roi));
            dataset(tr_n).ROI = ROI;
        end
        clear del_roi
        
        %committing motion registration parameters to disk to allow quicker
        %load times
        %regParams(dataset,'action','commit')       %don't use this line - it breaks ROI saving.
        

        %saving changes to the 2P object
        data = dataset;
        save([direc curr_set], 'data');
        clear data

        clear ROI
        
        
        
    else
        dataset = load([direc 'expt.mat']);
        dataset = dataset.data;
        curr_set = 'expt.mat';
    end
    
    
    %saving backup copy of raw data files before making any changes to raw
    %data files
    %backup_rawfiles(dataset);
    
    %use line below to recover rawfiles from backed-up copy
    %recover_rawfiles(dataset);
    
        
    %Attempting to automatically detect bad frames and bad trials and then blank 
    %them out. Bad frames are those that are dissimilar to the largest group of
    %similar frames
    %trial_quality_check(dataset, trial_n)
        
    %re-writing raw data file to reflect blanked out frames
    
    %write_to_rawfile(dataset, trial_n, new_matrix);
    
%---------------------------------------------
end
fclose(fid);


        
    
%-------------------------------------
fid = fopen(list_direc);
%loop to go through all experiment datasets listed in list file, once more to go through all steps that require manual intervention.
while 1
    
    direc = fgetl(fid);
       
    if ischar(direc) ~= 1
        break
    else
    end
    
    direc = [direc '\'];
    
    cd(direc)
    
    dataset = load([direc 'expt.mat']);
    dataset = dataset.data;
    n_trials = size(dataset, 2);
    
    
    
    
    %drawing individual cell ROIs 

    %identifying and getting rid of bad trials (usually trials with z motion in them)
    try 
        del = dataset(1).process.bad_trials;                    %checking if bad trials have already been identified
        %creating dataset_ROIs without bad trials so as to be able to draw
        %ROIs without being distracted by bad trials
        try
            good_trs = [];
            for trial_n = 1:n_trials
                
                if isempty(dataset(trial_n).stim) == 1
                    isgood = 1;
                else
                    isgood = dataset(trial_n).process.good_trial;
                end
                
                good_trs = [good_trs; isgood];
            end
            good_tr_list = find(good_trs == 0);
            
            dataset_ROIs = dataset(good_tr_list);
        catch
            keyboard
        end
    catch
        [good_trials, bad_trials] = find_bad_trials(dataset);
       
        dataset_ROIs = dataset(good_trials);                    %keeping only good trials for ROI-drawing function - will copy over ROIs to main dataset later
        
        %keeping track of which trials were manually identified as being
        %bad
        
        for trial_n = 1:n_trials
            if isempty(intersect(trial_n, good_trials)) == 1
                dataset(trial_n).process.good_trial = 1;
            else
            end
            if isempty(intersect(trial_n, bad_trials)) == 1
                dataset(trial_n).process.good_trial = 0;
            else
            end
        end
        
        dataset(1).process.bad_trials = 'done';
                
        %adding empty background ROIs to all trials
        ROI.level = 0.023;
        ROI.roi = zeros(255, 256) + 1;
        del_roi = abs(ROI.roi - 1);
        ROI.roi(1:10, 1:10) = 0;                  %tiny patch in top left corner to use as background for calculations
        ROI.notes = 'empty background ROI. DO NOT DELETE ME';
        for tr_n = 1:n_trials
            del = dataset(tr_n).imageStack;
            ROI.backgroundLevel = mean(mean(del(:, :, 10).*del_roi));
            dataset(tr_n).ROI = ROI;
        end
        clear del_roi
        
        
        %saving changes to the 2P object
        data = dataset;
        save([direc curr_set], 'data');
        clear data
        
        if (length(good_trials)./length(bad_trials) ) < 1
        disp(['Warning: ' int2str(round( (length(bad_trials)./n_trials).*100) ) ' percent of trials were bad.'])
        else
        end
        
    end
    
    %putting in option to manually choose to analyse pre-treatment and
    %post-treatment sets independently, i.e. not tracking the same cells
    %through both.
    try del = isempty(dataset(1).stim.switch_valve_trial);
        if isempty(dataset(1).stim.switch_valve_trial) == 0
            separate_sets = questdlg('Track same cells through baseline and treatment trials?',... 
                ' ',... 
                'Separate', 'Track through', 'Track through');

            switch separate_sets
                case 'Separate'
                    switch_trial = dataset(1).stim.switch_valve_trial;
                    dataset_baseline = dataset(1:switch_trial);
                    dataset_treatment = dataset((switch_trial + 1):size(dataset, 2) );
                case 'Track through'

            end

        else
            separate_sets = 'Track through';
        end
    catch
        separate_sets = 'Track through';
    end
    
    
    try
        switch separate_sets
            case 'Separate'
                dataset_baseline = selectROIs(dataset_baseline, 'soma');               %calling Rob's ROI-drawing function separately for baseline trials
                dataset(1:switch_trial) = dataset_baseline;
                dataset_treatment = selectROIs(dataset_treatment, 'soma');             %calling Rob's ROI-drawing function separately for treatment trials
                dataset( (switch_trial + 1):size(dataset, 2) ) = dataset_treatment;
            case 'Track through'
                dataset = selectROIs(dataset_ROIs, 'soma');               %calling Rob's ROI-drawing function
        end
    catch
        keyboard
    end
    
    cell_ROI_mat = dataset(1).ROI(2).roi;       %matrix of cell ROIs obtained from longest dataset
    
    
    %saving changes to the 2P object
    data = dataset;
    save([direc curr_set], 'data');
    clear data
    
    
    %checking cells manually for z-movt
    try
        bad_cells = ROI_z_movt_checker_20160505(dataset, direc);
    catch
        keyboard
    end
    
    %getting rid of bad cells from ROI mat
    ROI_mat_orig = dataset(1).ROI(2).roi;
    for bad_cell_n = 1:length(bad_cells);
        bad_cells = sort(bad_cells, 'descend');
        curr_bad_cell = bad_cells(bad_cell_n);
        curr_pixi = find(ROI_mat_orig == curr_bad_cell);
        ROI_mat_orig(curr_pixi) = 0;
        later_cellsi = find(ROI_mat_orig > curr_bad_cell);
        ROI_mat_orig(later_cellsi) = ROI_mat_orig(later_cellsi) -  1; 
    end
    
    for trial_n = 1:n_trials
        dataset(trial_n).ROI(2).roi = ROI_mat_orig;
    end
    
    
    %saving changes to the 2P object
    data = dataset;
    save([direc curr_set], 'data');
    clear data
    
    
end
fclose(fid);
    
%---------------------------------------------

fid = fopen(list_direc);
dir_counter = 1;
%loop to go through all experiment datasets listed in list file, once more to extract raw data traces from cell-ROIs.
while 1
    
    direc = fgetl(fid);
       
    if ischar(direc) ~= 1
        break
    else
    end
    
    direc = [direc '\'];
    
    cd(direc)
    
    dataset = load([direc 'expt.mat']);
    dataset = dataset.data;
    
    
    %checking if dataset was split for drawing ROIs or if cells were tracked through
    try
        switch_trial = dataset(1).stim.switch_valve_trial;
    catch
        switch_trial = [];
        
    end
    
    if isempty(switch_trial) == 0
        baseline_ROI = dataset(1).ROI(2).roi;
        treatment_ROI = dataset(switch_trial + 1).ROI(2).roi;
        if sum(sum(abs(baseline_ROI - treatment_ROI))) > 0
            split_set = 1;
        else
            split_set = 0;
        end
    else
        split_set = 0;
    end
    
    
    if split_set == 0
        %extracting raw, mean fluorescence traces for each cell-ROI and saving to a single file for further analysis.
        raw_data_mat = get_raw_traces(dataset);
        %saving extracted raw data
        save([direc 'expt_raw_traces.mat'], 'raw_data_mat');
    elseif split_set == 1       
        %if dataset is split for drawing ROIs, raw data is extracted separately for baseline and for treatment trials
        baseline_dataset = dataset(1:switch_trial);
        treatment_dataset = dataset( (switch_trial + 1):(size(dataset, 2)) );
        
        %extracting raw, mean fluorescence traces for each cell-ROI and saving to a single file for further analysis.
        raw_data_mat_baseline = get_raw_traces(baseline_dataset);
        %extracting raw, mean fluorescence traces for each cell-ROI and saving to a single file for further analysis.
        raw_data_mat_treatment = get_raw_traces(treatment_dataset);
        
        %saving extracted raw data
        mkdir([direc 'split_set'])
        save([direc 'split_set\expt_raw_traces_baseline.mat'], 'raw_data_mat_baseline');
        save([direc 'split_set\expt_raw_traces_treatment.mat'], 'raw_data_mat_treatment');
        
    else
    end
    
    
    disp(['Currently extracting data from dataset ' int2str(dir_counter)]);
    
    dir_counter = dir_counter + 1;
end
fclose(fid);

clear all
close all