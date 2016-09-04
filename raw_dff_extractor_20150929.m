clear all
close all

list_direc = ['D:\Data\CSHL\dataset_list_dopamine_20150421.txt'];
fid = fopen(list_direc);

%flip this switch to ignore 'dataset done' labels and re-analyse datasets
ignore_dones = 1;        %1 - ignore labels and re-analyse, 0 - don't re-analyse


direc_counter = 0;

%loop to go through all experiment datasets listed in list file
while 1
    direc_counter = direc_counter + 1;
    direc = fgetl(fid);
    %replacing C: with D:
    a = findstr(direc, 'C:');
    direc(a) = 'D';
    
       
    if ischar(direc) ~= 1
        break
    else
    end
    
    direc = [direc '\'];
    
    cd(direc)


    %importing raw .tiff files and saving 2P objects as .mat files; skips if
    %this has already been done before
    importSItree();
    

    %loading in 2P objects
    d=dir('*.mat');
    d(find(strcmp({d.name},'.')))=[];
    d(find(strcmp({d.name},'..')))=[];
    d(find(strcmp({d.name},'expt.mat')))=[];
    d(find(strcmp({d.name},'expt_raw_traces.mat')))=[];
    d(find(strcmp({d.name},'aved_tr_frames.mat')))=[];
    d(find(strcmp({d.name},'movie.mat')))=[]
    
    
    if exist([direc 'expt.mat']) ~= 2           %Note: expt is a reserved name, don't name the experiment folder expt.
        
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
    
    
%---------------------------------------------
end
fclose(fid);

clear dataset
        
    
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
    
    %replacing C: with D:
    a = findstr(direc, 'C:');
    direc(a) = 'D';
    
    
    cd(direc)
    
    
    dataset = load([direc 'expt.mat']);
    dataset = dataset.data;
    
    %skipping dataset if marked as manual processing 'done'.
    try del = dataset(1).process.manual_done;
        if del == 1
            if ignore_dones == 1        %not skipping and re-analysing because of manual over-ride
                continue
            else
            end
        else
        end
    catch
        
    end
    
    
    %identifying and getting rid of bad trials (usually trials with z motion in them)
    try 
        del = dataset(1).process.bad_trials;                    %checking if bad trials have already been identified
        if del == 'done';
            bad_tr_done = 1;
        else
            bad_tr_done = 0;
        end
    catch
        bad_tr_done = 0;
    end
       
        
        
    if bad_tr_done == 1
        %re-constructing dataset with only 'good trials' for ROI selection
        tr_counter = 1;
        good_trials = [];
        for trial_n = 1:size(dataset, 2);
            tr_type = dataset(trial_n).process.good_trial;
            if tr_type == 1
                dataset_ROIs(tr_counter) = dataset(trial_n);        %only good trials will be used to draw ROIs later.
                good_trials = [good_trials; trial_n];
                tr_counter = tr_counter + 1;
            else

            end

        end

        clear tr_counter
    else
         
        [good_trials, bad_trials] = find_bad_trials(dataset);
        dataset_ROIs = dataset(good_trials);                    %keeping only good trials for ROI-drawing function - will copy over ROIs to main dataset later
        n_trials = size(dataset, 2);
        
        %keeping track of which trials were manually identified as being
        %bad
        for trial_n = 1:n_trials
            if isempty(intersect(trial_n, bad_trials)) == 1
                dataset(trial_n).process.good_trial = 1;
            else
            end
            if isempty(intersect(trial_n, bad_trials)) == 0
                dataset(trial_n).process.good_trial = 0;
            else
            end
            
        end
        
       
        if (length(good_trials)./length(bad_trials) ) < 1
            disp(['Warning: ' int2str(round( (length(bad_trials)./n_trials).*100) ) ' percent of trials were bad.'])
        else
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
                    dataset_baseline = dataset_ROIs(1:switch_trial);
                    dataset_treatment = dataset_ROIs((switch_trial + 1):size(dataset_ROIs, 2) );
                case 'Track through'

            end

        else
            separate_sets = 'Track through';
        end
    catch
        separate_sets = 'Track through';
    end
    
    %try
        switch separate_sets
        case 'Separate'
            dataset_baseline = selectROIs_m(dataset_baseline, 'soma');               %calling Rob's ROI-drawing function separately for baseline trials
            dataset_ROIs(1:switch_trial) = dataset_baseline;
            dataset_treatment = selectROIs_m(dataset_treatment, 'soma');             %calling Rob's ROI-drawing function separately for treatment trials
            dataset_ROIs( (switch_trial + 1):size(dataset_ROIs, 2) ) = dataset_treatment;
        case 'Track through'
            dataset_ROIs = selectROIs_m(dataset_ROIs, 'soma');               %calling Rob's ROI-drawing function
        end
       
    
        
    %saving ROIs selected
    %transferring cell ROIs from dataset_ROIs (doesn't contain bad trials)
    %to dataset (contains tagged bad trials)
    ROI_mat = dataset_ROIs(1).ROI(2).roi;
    for trial_n = 1:size(dataset, 2)
        dataset(trial_n).ROI(2).roi = ROI_mat;
        dataset(trial_n).ROI(2).notes = 'soma';
        
    end
    
    %saving changes to the 2P object
    data = dataset;
    save([direc curr_set], 'data');
    clear data
        
    
    %looking at each cell's local neighbourhood across trials for movt
    [bad_cell_list, saved_lags] = bad_cell_catcher_20151023(dataset, direc, good_trials, 0);
    dataset(1).process.saved_lags = saved_lags;
    
    
    %removing bad cells from ROI matrix    
    ROI_mat = dataset_ROIs(1).ROI(2).roi;
    
    for bad_cell_n = 1:length(bad_cell_list)
        bad_cell_ni = bad_cell_list(bad_cell_n);
        
        del = find(ROI_mat == bad_cell_ni);
        ROI_mat(del) = 0;                   %getting rid of bad cell
        %subtracting 1 from the label for each subsequent ROI to keep
        %numbering continuous
        del = find(ROI_mat > bad_cell_ni);
        ROI_mat(del) = ROI_mat(del) - 1; 
    end
   
    
    %transferring cell ROIs from dataset_ROIs (doesn't contain bad trials)
    %to dataset (contains tagged bad trials)
    for trial_n = 1:size(dataset, 2)
        dataset(trial_n).ROI(2).roi = ROI_mat;
        dataset(trial_n).ROI(2).notes = 'soma';
        
    end
    
    
    %marking dataset as complete for manual processing
    all_done = input('All done? Stop re-opening dataset for manual processing? 1 - Yes, 0 - No');
    if all_done == 1
        dataset(1).process.manual_done = 1;
    else
    end
    
    
    %saving changes to the 2P object
    data = dataset;
    save([direc curr_set], 'data');
    clear data
    
   
end
fclose(fid);
    
%---------------------------------------------

clear dataset

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
    
    %replacing C: with D:
    a = findstr(direc, 'C:');
    direc(a) = 'D';
    
    cd(direc)
    
    dataset = load([direc 'expt.mat']);
    dataset = dataset.data;
    saved_lags = dataset(1).process.saved_lags;
    
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
    
    disp(['Currently extracting data from dataset ' int2str(dir_counter)]);
    if split_set == 0
        %extracting raw, mean fluorescence traces for each cell-ROI and saving to a single file for further analysis.
        
        raw_data_mat = get_raw_traces_20151026(dataset, direc, saved_lags);
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
    
    
    disp(['Extraction complete, moving to next set']);
    
    dir_counter = dir_counter + 1;
    
end
fclose(fid);
