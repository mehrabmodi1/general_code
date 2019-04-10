%this script visualises z-stacks acquired on the ScanImage rig.

clear all
close all

direc = 'D:\Data\CSHL\20160915\fly2\zstack\';
save_path = 'C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Talks\Labmeets\20161011\zstacks\';
saving = 1;

%reading in dataset parameters
try
    [n_frames, n_trials, frame_rate, zoom, d_path, f_name, tr_tg_num, ave_frame] = tiff_info(direc);
catch
    keyboard
end

n_slices = n_trials;


%reading in and averaging frames for each slice, and storing each averaged
%slice in a matrix for viewing as well as saving.
ave_mat = zeros(size(ave_frame, 1), size(ave_frame, 2), n_slices);
for slice_n = 1:n_slices
    for frame_n = 1:n_frames
        trial_no_f = slice_n + tr_tg_num - 1;           %adding the trial tag number of the first trial of this set to get the right filename
        trial_no_f = sprintf('%03.0f', trial_no_f);
        file_path = [direc f_name trial_no_f, '.tif'];  %full file path for currently analysed trial
        
        frame = double(imread(file_path, frame_n) );
        
        if frame_n == 1
            ave_frame = frame;
        else
        end
        
        if frame_n > 1
            ave_frame = ave_frame + frame;
        else
        end
        
    end
    ave_frame = ave_frame./n_frames;
    ave_mat(:, :, slice_n) = ave_frame;
end


%displaying and possibly saving averaged z-stack frames
for slice_n = 1:n_slices
    colormap('gray')
    figure(1)
    imagesc(ave_mat(:, :, slice_n))
    drawnow
    pause(0.1)
    
    %saving
    if saving == 1
        save([save_path 'direc.mat'], 'direc');
        print('-dtiff', [save_path 'slice_' int2str(slice_n) '.tif'])
        %imwrite(uint8(ave_mat(:, :, slice_n)), [save_path, 'zstack.tif'], 'WriteMode','append');
    else
    end
    
end

