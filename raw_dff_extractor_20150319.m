list_direc = ['C:\Data\CSHL\dataset_list_20150319.txt'];
fid = fopen(list_direc);

%loop to go through all experiment datasets listed in list file
while 1
    
    direc = fgetl(fid);
       
    if ischar(direc) ~= 1
        break
    else
    end
    
    direc = [direc '\'];
    
    %direc = 'C:\Data\CSHL\20150302\expt\';

    cd(direc)


    %importing raw .tiff files and saving 2P objects as .mat files; skips if
    %this has already been done before
    importSItree();

    %loading in 2P objects
    d=dir('*.mat');
    d(find(strcmp({d.name},'.')))=[];
    d(find(strcmp({d.name},'..')))=[];


    %identifying longest dataset in curr expt folder to use as basis of identifying individual cell ROIs
    t_vec = zeros(1, length(d));
    for set_n = 1:length(d)
        curr_set = d(set_n).name;

        %loading 2P object
        dataset = load(curr_set);
        dataset = dataset.data;
        n_trials = size(dataset, 2);
        t_vec(set_n) = n_trials;
    end
    [del long_set] = max(t_vec);
    clear del
    clear t_vec

    %loop to prep all datasets for data extraction before using the longest
    %dataset to indentify single cell ROIs
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

            %aligning the aligned trials to each other
            alignRepeats(dataset)
            
            
            
            %updating info in 2P object
            updateMuStack(dataset)
            
            
            %saving changes to the 2P object
            data = dataset;
            save([direc curr_set], 'data');
            clear data
            
         end
       
       
        
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
        
        clear ROI

     cd(direc)

     
    end    
    
    %drawing individual cell ROIs using the single, longest dataset in curr
    %expt folder
    curr_set = d(long_set).name;

    %loading 2P object
    dataset = load(curr_set);
    dataset = dataset.data;
    n_trials = size(dataset, 2);

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
        
    dataset = selectROIs(dataset, 'soma')               %calling Rob's ROI-drawing function
    disp('Long set ROIs')
    cell_ROI_mat = dataset(1).ROI(2).roi;       %matrix of cell ROIs obtained from longest dataset
    
    %saving changes to the 2P object
    data = dataset;
    save([direc curr_set], 'data');
    clear data
    
    
    %looping through all datasets in current expt folder and applying the
    %same cell-ROI matrix
    for set_n = 1:length(d)
        curr_set = d(set_n).name;

        %loading 2P object
        dataset = load(curr_set);
        dataset = dataset.data;
        n_trials = size(dataset, 2);
        
        %copying over ROIs drawn for longest trial set in this experiment (see above)
        for trial_n = 1:n_trials
            dataset(trial_n).ROI(2).roi = cell_ROI_mat;
            dataset(trial_n).ROI(2).notes = 'soma';
        end
        
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
        
        
        %running Rob's ROI selection algorithm to re-check ROIs
        dataset = selectROIs(dataset, 'soma')       %calling Rob's ROI-drawing function
        disp(['ROIs for set ' int2str(set_n)])
        
        %saving changes to the 2P object
        data = dataset;
        save([direc curr_set], 'data');
        clear data
               
        raw_data_mat = get_raw_traces(dataset);
        keyboard
    end


end

fclose(fid);

