function [raw_data_mat] = get_raw_traces(dataset)

%extracting raw intensity data, averaged over ROI pixels and saving
%to file for further analysis
cell_roi_mat = dataset(1).ROI(2).roi;

%identifying centroids of each ROI
centroids = regionprops(cell_roi_mat, 'Centroid');
centroids = round(cat(1, centroids.Centroid));

n_cells = max(max(cell_roi_mat));
n_trials = size(dataset, 2);

%allowing for different numbers of frames in each trial by setting n_frames
%to the largest number of frames of any trial in the dataset.
n_frames_list = zeros(1, n_trials);
bad_tr_list = [];
for trial_n = 1:n_trials
    if isempty(dataset(trial_n).stim) == 1
        bad_tr_list = [bad_tr_list; trial_n]; 
    else
        n_frames_list(1, trial_n) = size(dataset(trial_n).info.micronsPerPixel_XAxis, 2);
    end
end
good_tr_list = 1:n_trials;
good_tr_list(bad_tr_list) = [];
dataset = dataset(good_tr_list);
n_frames_list(bad_tr_list) = [];
n_trials = length(good_tr_list);

n_frames_max = max(n_frames_list);
raw_data_mat = zeros(n_frames_max, n_cells, n_trials) + nan;

for trial_n = 1:n_trials

    n_frames = n_frames_list(1, trial_n);           %number of frames may differ from trial to trial
    curr_stack = dataset(trial_n).imageStack;
    
        
    for cell_n = 1:n_cells
        curr_roi = zeros(size(cell_roi_mat, 1), size(cell_roi_mat, 2) );
        del = find(cell_roi_mat == cell_n);
        curr_roi(del) = 1;                             %ROI of current cell
        for frame_n = 1:n_frames
            raw_data_mat(frame_n, cell_n, trial_n) = nanmean(nanmean(curr_stack(:, :, frame_n).*curr_roi) );
            
        end
    end
end