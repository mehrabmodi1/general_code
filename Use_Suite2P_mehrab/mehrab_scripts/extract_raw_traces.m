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

%checking if any trials have already been extracted and recovering
if exist([save_path, 'extracted_raw_data_mat.mat']) == 2
    try
        raw_data_mat = load([save_path, 'extracted_raw_data_mat.mat']);
        raw_data_mat = raw_data_mat.raw_data_mat;
        done_trs = find(squeeze(isnan(raw_data_mat(1, 1, :))) == 0);
        start_trial = done_trs(end) + 1;
        disp(['Recovered some extracted traces. Starting to extract trial ', int2str(start_trial), '.'])

        ref_im = load([save_path, 'ref_im.mat']);
        ref_im = ref_im.ref_im;

        time_stamps = load([save_path, 'tif_time_stamps.mat']);
        time_stamps = time_stamps.time_stamps;
    catch
        start_trial = 1;
    end
else
    start_trial = 1;
end

for trial_n = start_trial:n_trials
    try
        %reading in stack object
        im_obj = ScanImageTiffReader([direc, dir_contents(trial_n).name]);
        %obtaining image stack
        stack = im_obj.data();
        stack = permute(stack,[2 1 3]);
        %obtaining, logging timestamp
        curr_time_stamp = parse_tiff_timestamp(im_obj);
        time_stamps(trial_n).tstamp = curr_time_stamp;
        time_stamps(trial_n).name = dir_contents(trial_n).name;
    catch
        continue
    end
    
    %applying a slow x-y motion correction assuming negligible motion
    %within a trial
    if start_trial == 1 && trial_n == 1
        ref_im = mean(stack, 3, 'omitnan');
    else
        stack = slow_xy_motion_correct(stack, ref_im);
    end

    %reading in raw fluorescence data for each ROI into a single matrix (dim1 - frame_n, dim2 - ROI_n)
    raw_vec = zeros(1, n_cells);
    n_frames = size(stack, 3);
    
    for frame_n = 1:n_frames
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
    
    %saving extracted traces
    if exist(save_path) == 0
        mkdir(save_path)
    else
    end
    save([save_path, 'extracted_raw_data_mat.mat'], 'raw_data_mat');
    save([save_path, 'ref_im.mat'], 'ref_im');
    save([save_path, 'tif_time_stamps.mat'], 'time_stamps');
    disp(['traces extracted, from trial ', int2str(trial_n), ', and saved.'])

end
cd(prev_direc)