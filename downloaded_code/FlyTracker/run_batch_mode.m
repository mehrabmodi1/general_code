clear all
close all

%Note: Parameters need to be set manually, folders in list should run one
%by one, assuming calibration has already been done.

curr_batch_folder = 'E:\Data\Raw_Data_Current\behavior_movies\2021-11-03\';     %folder that contains one folder for each experiment (ie. a single calibration run)
expt_list = dir(curr_batch_folder);
expt_list(1:2) = [];

folders = {''};
%looping through folders for manual calibration of each experiment.
for dir_n = 1:size(expt_list, 1)
    curr_dir = [curr_batch_folder, expt_list(dir_n).name, '\'];      %path for current experiment
    folders{dir_n} = curr_dir;
    %checking if curr_expt is already clibrated
    if exist([curr_dir, 'calibration.mat']) == 2
        continue
    else
        cd(curr_dir);
        disp(curr_dir);
        disp('Open current folder and clibrate manually, then close tracker and enter ''dbcont''. ');
        tracker;        
        keyboard
    end
end


%running tracking

% set options (omit any or all to use default options)
%options.granularity  = 10000;
options.num_chunks   = 16;       % set either granularity or num_chunks
options.num_cores    = 8;
options.max_minutes  = Inf;
options.save_JAABA   = 1;
options.save_seg     = 1;

% loop through all folders
for f=1:numel(folders)
    % set parameters for specific folder
    videos.dir_in  = folders{f};
    videos.dir_out = folders{f}; % save results into video folder
    videos.filter = '*.ufmf';     % extension of the videos to process
    
    % track all videos within folder
    tracker(videos,options);
end