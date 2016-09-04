direc = 'C:\Data\CSHL\20150302\expt\';

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
        alignStack(dataset, [1,20], 'verbose', 1)

        %updating info in 2P object
        updateMuStack(dataset)

        %aligning the aligned trials to each other
        try
            alignRepeats(dataset)
        catch
            keyboard
        end
        
        %saving changes to the 2P object
        data = dataset;
        save([direc curr_set], 'data');
        clear data
    end
    
    %identification of background ROI
    if isempty(dataset(1).ROI) == 1
        dataset = ROI_batch(dataset)
        %saving changes to the 2P object
        data = dataset;
        save([direc curr_set], 'data');
        clear data
    else
        disp('background ROI already added.')
    end
    
    

 cd(direc)
    
end    
   


%drawing individual cell ROIs using the single, longest dataset in curr
%expt folder
curr_set = d(long_set).name;

%loading 2P object
dataset = load(curr_set);
dataset = dataset.data;
n_trials = size(dataset, 2);

%identifying and getting rid of bad trials (usually trials with z motion in
%them
[im, del, im_stack] = visualiseDrift(dataset, [], 'mean');
c=corrcoef(im);                 %corrcoef matrix of mean baseline frame of each trial, across trials
del = eye(size(c, 1) );         
del = find(del == 1);           %getting rid of diagonal 1's
c(del) = nan;
clear del

c = nanmean(c);                 %vector of corrcoefs of each trial with all other trials
corr_cutoff = c(1);             %cutoff to get rid of non-matching trials is the mean corrcoef of tr1 with all other trials

bad_trials = find(c < corr_cutoff);
good_trials = 1:n_trials;
good_trials(bad_trials) = [];   %list of trials with mean corrcoef higher than cutoff defined above

dataseta = dataset(good_trials);

keyboard
dataset = selectROIs(dataset(1))       %calling Rob's ROI-drawing function

%saving changes to the 2P object
data = dataset;
save([direc curr_set], 'data');
clear data

cell_ROI_mat = dataset(1).ROI(2).roi;       %matrix of cell ROIs obtained from longest dataset


%looping through other datasets in current expt folder and applying the
%same cell-ROI matrix
for set_n = 1:length(d)
    curr_set = d(set_n).name;

    %loading 2P object
    dataset = load(curr_set);
    dataset = dataset.data;
    n_trials = size(dataset, 2);

    for trial_n = 1:n_trials
        dataset(1).ROI(2).roi = cell_ROI_mat;
    end

    %running Rob's ROI selection algorithm to re-check ROIs
    dataset = selectROIs(dataset)       %calling Rob's ROI-drawing function
    
    keyboard
    %saving changes to the 2P object
    data = dataset;
    save([direc curr_set], 'data');
    clear data
end

% for tr = 1:40
%     figure(1)
%     imagesc(im(:, :, tr));
%     colormap('gray');
%     title(int2str(tr));
%     a = input('press enter')
% end
