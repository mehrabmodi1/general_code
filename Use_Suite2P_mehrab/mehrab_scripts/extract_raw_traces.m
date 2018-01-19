function [raw_data_mat] = extract_raw_traces(direc, ROI_mat, save_path)
%syntax: [raw_data_mat] = extract_raw_traces(direc, ROI_mat, save_path)
%This function extracts raw traces from all the large tif files in direc using 
%the specified ROI_mat. It saves the extracted raw_data_mat after each trial 
%is extracted and picks up where it left off.

prev_direc = pwd;
cd(direc)
dir_contents = dir('*.tif');
n_trials = size(dir_contents, 1);
remove_small_tifs([direc]);
n_cells = size(ROI_mat, 3);
raw_data_mat = zeros(100, n_cells, n_trials) + nan;
n_frames = zeros(n_trials, 1);

for trial_n = 1:n_trials
    stack = ScanImageTiffReader([direc, dir_contents(trial_n).name]).data();
    stack = permute(stack,[2 1 3]);
    %applying a slow x-y motion correction assuming negligible motion
    %within a trial
    if trial_n == 1
        ref_im = mean(stack, 3, 'omitnan');
    else
        stack = slow_xy_motion_correct(stack, ref_im);
    end

    %reading in raw fluorescence data for each ROI into a single matrix (dim1 - frame_n, dim2 - ROI_n)
    raw_vec = zeros(1, n_cells);
    n_frames = size(stack, 3);
    
    for frame_n = 1:(100 + round(rand(1, 1).*10))
        curr_frame = stack(:, :, frame_n);
        curr_frame = im2double(curr_frame);

        for ROI_n = 1:n_cells
            curr_ROI = ROI_mat(:, :, ROI_n);
            %curr_ROI_pix = find(curr_ROI == 1);
            raw_vec(1, ROI_n) = mean(curr_frame(curr_ROI == 1));

        end
        curr_raw_data_mat(frame_n, :) = raw_vec;   %raw data mat for current trial
        
    end
    
    %putting raw data mat for current trial into big matrix for all trials
    diff_n_frames = size(raw_data_mat, 1) - size(curr_raw_data_mat, 1);
    if sign(diff_n_frames) == 1
        raw_data_mat(1:size(curr_raw_data_mat, 1), :, trial_n) = curr_raw_data_mat;
    elseif sign(diff_n_frames) == -1
        raw_data_mat_orig = raw_data_mat;
        raw_data_mat = zeros(size(curr_raw_data_mat, 1), n_cells, n_trials) + nan;
        raw_data_mat(1:size(raw_data_mat_orig, 1), :, :) = raw_data_mat_orig;
        raw_data_mat(1:size(curr_raw_data_mat, 1), :, trial_n) = curr_raw_data_mat;
       
    else
    end
    
    clear stack
    
    %PICK UP THREAD HERE
    %implement saving extracted traces and recovery in the case of
    %interruptions.
    disp(['traces extracted, from trial ', int2str(trial_n), ', and saved.'])
    
end

keyboard
cd(prev_direc)