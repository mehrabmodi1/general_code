function [raw_data_mat] = extract_raw_traces(direc, ROI_mat, save_path, do_registration, extract_both_channels, subtraction_spec_vec)
%syntax: [raw_data_mat] = extract_raw_traces(direc, ROI_mat, save_path)
%This function extracts raw traces from all the large tif files in direc using 
%the specified ROI_mat. It saves the extracted raw_data_mat after each trial 
%is extracted and picks up where it left off.

prev_direc = pwd;
cd(direc)
try
    dir_contents = dir_date_sorted(direc, '*.tif');
catch
    keyboard
end
n_trials = length(dir_contents);
remove_small_tifs(direc);
n_cells = size(ROI_mat, 3);
raw_data_mat = zeros(100, n_cells, n_trials) + nan;
n_frames = zeros(n_trials, 1);



%checking if any trials have already been extracted and recovering
if exist([save_path, 'extracted_raw_data_mat.mat']) == 2
    try
        raw_data_mat = load([save_path, 'extracted_raw_data_mat.mat']);
        raw_data_mat = raw_data_mat.raw_data_mat;
        
        if length(size(raw_data_mat)) == 3
            done_trs = find(squeeze(isnan(raw_data_mat(1, 1, :))) == 0);
            if isempty(done_trs) == 1
                done_trs = 0;
            else
            end
        else
            done_trs = 0;
        end
        
        if max(done_trs) == n_trials
            disp('all traces already extracted... skipping.')
            del = [];
            save([save_path, 'trace_extraction_complete.mat'], 'del');
           
            return
        else
            start_trial = max(done_trs) + 1;
           
        end
        
        
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

%reading in previously saved background region ROI for background (PMT offset) subtraction
if exist([save_path, '\bk_ROI.mat']) == 2
    bk_ROI = load([save_path, '\bk_ROI.mat']);
    bk_ROI = bk_ROI.bk_ROI;
else
    bk_ROI = [];
    keyboard
end


for trial_n = start_trial:n_trials
    try
        %reading in stack object
        im_obj = ScanImageTiffReader([direc, dir_contents(trial_n).name]);
        %obtaining image stack
        stack = im_obj.data();
        %stack_orig = stack;
        stack = permute(stack,[2 1 3]);
        
        %obtaining, logging timestamp
        curr_time_stamp = parse_tiff_timestamp(im_obj);
        time_stamps(trial_n).tstamp = curr_time_stamp;
        time_stamps(trial_n).name = dir_contents(trial_n).name;
        
        %checking how many color channels there are
        [frame_time, zoom, n_chans] = SI_tif_info(im_obj);
        if extract_both_channels == 1
            n_chans = 1;
        else
        end
        
        if n_chans == 2
            stack = stack(:, :, 1:2:end);       %getting rid of red channel frames
        else
        end
  
        clear frame_time
        clear zoom
    catch
        keyboard
        disp('skipping trace extraction from unreadable trials...')
        continue            %skipping trace extraction from unreadable trials - this doesn't matter because time-stamps will be matched up later with the stim param file anyway.
    end
    
%     p = gcp('nocreate');
%     if isempty(p) == 1
%         parpool(7)
%     else
%     end
    
    %applying a slow x-y motion correction assuming negligible motion
    %within a trial, unless otherwise specified by user
    if start_trial == 1 && trial_n == 1
        ref_im = mean(stack, 3, 'omitnan');
        %ROI_mat_curr = ROI_mat;
        if do_registration == 1
            saved_lags = zeros(size(stack, 3), 2);
            stack_reg = stack;
        else
        end
        %doing per-frame motion correction if specified by user
        if do_registration == 2
            %checking if detailed lags exist for current trial from a previous analysis run
            if exist([save_path, 'detailed_xy_lags.mat']) == 2
                detailed_lags = load([save_path, 'detailed_xy_lags.mat']);
                detailed_lags = detailed_lags.detailed_lags;
                local_registration_type = 3;
                detailed_lag_mat = detailed_lags(trial_n).frame_lags;
            
            else
                local_registration_type = do_registration;
                detailed_lag_mat = [];
            end

            try
                
                disp(['registering frames for trial ' int2str(trial_n)])
                [stack_reg, saved_lags] = slow_xy_motion_correct(stack, ref_im, [0, 0], detailed_lag_mat, local_registration_type);        %lag_mat is [0, 0] because this is trial 1.
                disp(['registration done. extracting data.'])
            catch
                keyboard
            end
            
        end
    else
        
        if do_registration == 1 %assume negligible motion within a trial and correct all frames within a stack with the same, manually determined lags
           	
            if exist([save_path, '\xy_lags.mat']) == 2
                lag_mat = load([save_path, '\xy_lags.mat']);
                
                lag_mat = lag_mat.lag_mat;
                
            else
                lag_mat = [];
            end
            
            %checking if detailed lags exist for current trial from a previous analysis run
            if exist([save_path, 'detailed_xy_lags.mat']) == 1
                detailed_lags = load([save_path, 'detailed_xy_lags.mat']);
                detailed_lags = detailed_lags.detailed_lags;
                local_registration_type = 3;
                detailed_lag_mat = detailed_lags(trial_n).frame_lags;
            
            else
                local_registration_type = do_registration;
                detailed_lag_mat = [];
            end
            
            disp(['registering frames for trial ' int2str(trial_n)])
            [stack_reg, saved_lags] = slow_xy_motion_correct(stack, ref_im, lag_mat(trial_n, :), detailed_lag_mat, local_registration_type);
            disp(['registration done. extracting data.'])
            
            
        elseif do_registration == 2 %do a 2D cc to identify lags for each frame in a stack relative to the meam frame for that stack after doing the slow correction above
            if exist([save_path, '\xy_lags.mat']) == 2
                lag_mat = load([save_path, '\xy_lags.mat']);
                
                lag_mat = lag_mat.lag_mat;
                
            else
                lag_mat = [];
            end
            
            %checking if detailed lags exist for current trial from a previous analysis run
            if exist([save_path, 'detailed_xy_lags.mat']) == 2
                detailed_lags = load([save_path, 'detailed_xy_lags.mat']);
                detailed_lags = detailed_lags.detailed_lags;
                local_registration_type = 3;
                detailed_lag_mat = detailed_lags(trial_n).frame_lags;
            
            else
                local_registration_type = do_registration;
                detailed_lag_mat = [];
            end
            
            disp(['registering frames for trial ' int2str(trial_n)])
            try
                [stack_reg, saved_lags] = slow_xy_motion_correct(stack, ref_im, lag_mat(trial_n, :), detailed_lag_mat, local_registration_type);
            catch
                keyboard
            end
            disp(['registration done. extracting data.'])
        end
                
    end
    
    detailed_lags(trial_n).frame_lags = saved_lags;
   
    %% reading in raw fluorescence data for each ROI into a single matrix (dim1 - frame_n, dim2 - ROI_n)
    n_frames = size(stack, 3);
   
    %checking for line noise and subtracting it away.
    if isa(stack_reg, 'int16') == 1
        stack_reg = double(stack_reg);
    else
    end   
        
    for frame_n = 1:n_frames
        curr_fr = stack_reg(:, :, frame_n);
        if subtraction_spec_vec(1, 2) == 1
            bk_fr = check_sig_noise(curr_fr);
        elseif subtraction_spec_vec(1, 2) == 0
            bk_fr = zeros(size(curr_fr, 1), size(curr_fr, 2));
        else
        end
        
        if median(median(bk_fr)) < median(median(curr_fr)).*0.05
            curr_fr_a = curr_fr - bk_fr;
            stack_reg(:, :, frame_n) = curr_fr_a;
        else
        end

    end
   
    %registration testing lines of code
%     if trial_n > 1 && rem(trial_n, 5) == 0
%         figure(2)
%         subplot(2, 1, 1)
%         imagesc(std(stack_reg, [], 3, 'omitnan'), [0, 10])
%         colormap('gray')
%         ROI_im = sum(ROI_mat, 3);
%         hold on
%         s = imagesc(ROI_im);
%         alpha(s, 0.5);
%         
%         subplot(2, 1, 2)
%         imagesc(std(stack_reg, [], 3, 'omitnan'))
%         colormap('gray')
%         
%         figure(3)
%         try
%             diff_fr = std(stack_reg, [], 3, 'omitnan') - first_fr;
%             imagesc(diff_fr, [-10, 10]);
%             colormap('gray')
%         catch
%             first_fr = std(stack_reg, [], 3, 'omitnan');
%             diff_fr = std(stack_reg, [], 3, 'omitnan') - first_fr;
%             imagesc(diff_fr, [-10, 10]);
%             colormap('gray')
%         end
%     keyboard
%     elseif trial_n == 1
%         first_fr = std(stack_reg, [], 3, 'omitnan');
%     end
%     
    
    bk_pixi = find(bk_ROI == 1);
    
    if isempty(bk_pixi)     %bk_ROI not saved, re-drawing
        disp('need to re-draw bk_ROI')
        imagesc(stack_reg(:, :, 1))
        bk_ROI = roipoly();
        save([save_path, '\bk_ROI.mat'], 'bk_ROI');
        bk_pixi = find(bk_ROI == 1);
        
    else
    end
    
    if subtraction_spec_vec(1, 1) == 1
        %computing a moving-window averaged bk_val for each frame.
        if isempty(bk_pixi) == 0

            bk_val_vec = zeros(size(stack_reg, 3), 1) + nan;
            for frame_n = 1:size(stack_reg, 3)
                curr_frame = stack_reg(:, :, frame_n);
                curr_frame = double(curr_frame);

                bk_val_vec(frame_n, 1) = mean(mean(curr_frame(bk_pixi)));

            end
            %computing moving-window average of bk_vals
            bk_val_vec = movmean(bk_val_vec, 101);

        else
        end
    else
    end
    
    for frame_n = 1:size(stack_reg, 3)
        
        curr_frame = stack_reg(:, :, frame_n);
        curr_frame = double(curr_frame);
        
        %doing a background subtraction to subtract away the PMT offset, if any
        if subtraction_spec_vec == 1
            if isempty(bk_pixi) == 1
                bk_val = 0;
            else
                bk_val = bk_val_vec(frame_n, 1);
            end

            curr_frame = curr_frame - bk_val;
        else
        end
        raw_vec = zeros(1, n_cells);
        
        for ROI_n = 1:n_cells
            try
                curr_ROI = ROI_mat(:, :, ROI_n);
            catch
                keyboard
            end
                    
            raw_vec(1, ROI_n) = mean(curr_frame(curr_ROI == 1));
            %disp(['extracting data from ROI ' int2str(ROI_n)]); 
        end
       
        curr_raw_data_mat(frame_n, :) = raw_vec;   %raw data mat for current trial
        %disp(['extracting data from frame ' int2str(frame_n)]); 
    end
    
    %putting raw data mat for current trial into big matrix for all trials
    diff_n_frames = size(raw_data_mat, 1) - size(curr_raw_data_mat, 1);
    if sign(diff_n_frames) == 1
        raw_data_mat(1:size(curr_raw_data_mat, 1), :, trial_n) = curr_raw_data_mat;
    elseif sign(diff_n_frames) == -1
        raw_data_mat_orig = raw_data_mat;
        raw_data_mat = zeros(size(curr_raw_data_mat, 1), n_cells, n_trials) + nan;
        try
            
            raw_data_mat(1:size(raw_data_mat_orig, 1), :, :) = raw_data_mat_orig;
        catch
            keyboard
        end
        raw_data_mat(1:size(curr_raw_data_mat, 1), :, trial_n) = curr_raw_data_mat;
       
    else
        raw_data_mat(:, :, trial_n) = curr_raw_data_mat;
    end
    
    clear stack
    clear stack_reg
    clear curr_raw_data_mat
    
    %saving extracted traces
    if exist(save_path) == 0
        mkdir(save_path)
    else
    end
    save([save_path, 'extracted_raw_data_mat.mat'], 'raw_data_mat');
    save([save_path, 'ref_im.mat'], 'ref_im');
    save([save_path, 'tif_time_stamps.mat'], 'time_stamps');
    disp(['traces extracted, from trial ', int2str(trial_n), ', and saved.'])
   
    try
        save([save_path, '\detailed_xy_lags.mat'], 'detailed_lags')
    catch
        keyboard
    end
        
end



del = [];
save([save_path, 'trace_extraction_complete.mat'], 'del');

cd(prev_direc)