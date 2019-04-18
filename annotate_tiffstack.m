clear all
close all

direc = 'C:\Users\Mehrab\Google Drive\Backup\Stuff\Janelia\Talks\DIG\20171019\sus_act_movies\20180626_fly3\';
fname = 'odor_trs_00009.tif';

stk = ScanImageTiffReader([direc, fname]).data();
stk = double(stk);
max_val = max(max(max(stk)));
%stk = stk./max_val;   %normlising frame values
moving_ave_window = 20;
stim_frames = [257, 857] - moving_ave_window;       %pairs of onset frame, off frame

tif_size = [size(stk, 1), size(stk, 2)];
n_frames = size(stk, 3);

stim_sq_ROI = zeros(tif_size(1), tif_size(2));
stim_sq_ROI( (tif_size(1).*0.8):(tif_size(1).*0.9), (tif_size(2).*0.8):(tif_size(2).*0.9)) = 1;
stim_sq_ROI_i = find(stim_sq_ROI == 1);

n_stims = size(stim_frames, 1);
ave_fr_mat = zeros(tif_size(1), tif_size(2), moving_ave_window);         %emtpy matrix to keep the moving window of frames in memory
averaged_stack = zeros(tif_size(1), tif_size(2), n_frames);              %empty matrix to start storing moving window averaged frames when they become available
first_frame = 1;

%creating mean baseline image for difference image stack
baseline_im = mean(stk(:, :, 1:stim_frames(1, 1)), 3);

for frame_n = 1:n_frames
    curr_frame = stk(:, :, frame_n);
    curr_frame_diff = curr_frame - baseline_im;
    
    for stim_n = 1:n_stims
        curr_stim_frames = stim_frames(stim_n, :);
    
        %checking if current fr is a stimulus fr and adding stimulus marker
        if frame_n >= curr_stim_frames(1) && frame_n <= curr_stim_frames(2)
            curr_frame(stim_sq_ROI_i) = max_val.*1.2;                   %stim indicator is 20% brighter than brightest pixel in stack
            curr_frame_diff(stim_sq_ROI_i) = max_val.*1.2;
        else
        end
    end
    
        %loop to implement moving window averaging of frames
    if frame_n <= moving_ave_window
        ave_fr_mat(:, :, frame_n) = curr_frame;
        ave_fr_mat_diff(:, :, frame_n) = curr_frame_diff;
        
    elseif frame_n > moving_ave_window
        ave_fr = mean(ave_fr_mat, 3, 'omitnan');
        averaged_stack(:, :, (frame_n - moving_ave_window)) = ave_fr;
        ave_fr_diff = mean(ave_fr_mat_diff, 3, 'omitnan');
        averaged_stack_diff(:, :, (frame_n - moving_ave_window)) = ave_fr_diff;
        
        ave_fr_mat(:, :, 1) = [];
        ave_fr_mat = cat(3, ave_fr_mat, curr_frame);
        ave_fr_mat_diff(:, :, 1) = [];
        ave_fr_mat_diff = cat(3, ave_fr_mat_diff, curr_frame_diff);
        
        %constructing and saving a combined frame for raw and diff frames
        width = size(ave_fr, 2);
        height = size(ave_fr, 1);
        col_width = round(height.*0.05);
        ave_fr_combined = zeros((2.*height + col_width), width) + max_val.*1.2;
        ave_fr_combined(1:width, 1:height) = ave_fr;
        ave_fr_combined((width + col_width + 1):size(ave_fr_combined, 1), 1:height) = ave_fr_diff;
        
        if first_frame == 1
            imwrite(uint8(ave_fr)', [direc, 'annotated_stk.tif']);
            imwrite(uint8(ave_fr_diff)', [direc, 'annotated_diff_stk.tif']);
            imwrite(uint8(ave_fr_combined)', [direc, 'annotated_combined_stk.tif']);
            first_frame = 0;
        elseif first_frame == 0
            imwrite(uint8(ave_fr)', [direc, 'annotated_stk.tif'], 'WriteMode', 'append');
            imwrite(uint8(ave_fr_diff)', [direc, 'annotated_diff_stk.tif'], 'WriteMode', 'append');
            imwrite(uint8(ave_fr_combined)', [direc, 'annotated_combined_stk.tif'], 'WriteMode', 'append');
        end
        
        
    else
    end
        
        
disp(frame_n)
    
end