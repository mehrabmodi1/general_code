function [] = tiff_downsampler(direc, t_factor)
%Syntax: function [] = tiff_downsampler(direc, t_factor)
%this function opens all the .tiff files in direc and downsamples them by
%averaging every t_factor frames and saving a downsampled copy. The
%original file with all the frames is then deleted.

% direc = 'D:\Data\Suite2P_registered\20170313\MB149BX20xsytGC6f - Copy\odor_trials\Plane1';
%t_factor = 10;

prev_direc = pwd;
cd(direc);
dir_contents = dir([direc '\*.tif']);
n_trials = size(dir_contents, 1);


for trial_n = 1:n_trials
    %loading in frames
    fname = dir_contents(trial_n).name;
    
    %skipping current file if it has already been down-sampled
    if isempty(findstr(fname, 'd_sampled')) == 0
        continue
    else
    end
    
    file_path = [direc '\' fname];
    stack = ScanImageTiffReader(file_path).data();
    n_frames = size(stack, 3);
    
    %avraging every t_factor frames and creating a stack sub-sampled in time
    n_ave_frames = ceil(n_frames./t_factor) - 1;        %number of down-sampled frames
    avestack = zeros(size(stack, 1), size(stack, 2), n_ave_frames) + nan;
    for ave_frame_n = 0:(n_ave_frames)
        if ave_frame_n < n_ave_frames
            curr_frames = stack(:, :, ( (ave_frame_n.*t_factor) + 1):( (ave_frame_n + 1).*t_factor));
        elseif ave_frame_n == n_ave_frames
            curr_frames = stack(:, :, ( (ave_frame_n.*t_factor) + 1):(size(stack, 3)));
        end
        ave_frame = nanmean(curr_frames, 3);
        avestack(:, :, (ave_frame_n + 1)) = ave_frame;
        ave_frame = uint16(ave_frame);
        
        %saving down-sampled tiff file to disk
        if ave_frame_n == 0
            imwrite(ave_frame, [file_path(1:(end-4)), '_d_sampled.tif']);
            
        elseif ave_frame_n > 0
            imwrite(ave_frame, [file_path(1:(end-4)), '_d_sampled.tif'], 'WriteMode', 'append');
        end
        
    end
    
    delete(file_path);
    
           
end


cd(prev_direc);

end