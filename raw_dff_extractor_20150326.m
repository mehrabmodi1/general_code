clear all
close all

list_direc = ['C:\Data\dataset_list_glc_ara_20150409.txt'];
fid = fopen(list_direc);

%loop to go through all experiment datasets listed in list file
while 1
    
    direc = fgetl(fid);
       
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
        regParams(dataset,'action','commit')

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


    %drawing individual cell ROIs 

    %identifying and getting rid of bad trials (usually trials with z motion in them)
    try 
        del = dataset(1).process.bad_trials;                    %checking if bad trials have already been identified
    catch
        [good_trials, bad_trials] = find_bad_trials(dataset);
        dataset = dataset(good_trials);
        dataset(1).process.bad_trials = 'done';
        
        %saving changes to the 2P object
        data = dataset;
        save([direc curr_set], 'data');
        clear data
        
        if (length(good_trials)./length(bad_trials) ) < 1
        disp(['Warning: ' int2str(round( (length(bad_trials)./n_trials).*100) ) ' percent of trials were bad.'])
        else
        end
        
    end
    
    try
        dataset = selectROIs(dataset, 'soma')               %calling Rob's ROI-drawing function
    catch
        keyboard
    end
    
    cell_ROI_mat = dataset(1).ROI(2).roi;       %matrix of cell ROIs obtained from longest dataset
    
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

    %extracting raw, mean fluorescence traces for each cell-ROI and saving to a single file for further analysis.
    raw_data_mat = get_raw_traces(dataset);
    %saving extracted raw data
    save([direc 'expt_raw_traces.mat'], 'raw_data_mat');
    
    
    
    disp(['Currently extracting data from dataset ' int2str(dir_counter)]);
    
    dir_counter = dir_counter + 1;
end
fclose(fid);

clear all
close all