function [stack] = load_all_frames(fpath)
%this function uses the matlab imread function in a loop to read in the
%all frames in a tif image. there is no need to specify the number of
%frames to read before-hand.

%path = 'E:\Data\Mehrab\20170328\MB594_opGC6f_35percent_laser\1\25_percent_power_00001.tif';

frame = imread(fpath, 1);
stack = zeros(size(frame, 1), size(frame, 2), 500) + nan;
stack(:, :, 1) = frame;
frame_n = 1;
while 1
    frame_n = frame_n + 1;
    if rem(frame_n, 500) == 0
        blank_frames = zeros(size(stack, 1), size(stack, 2), (size(stack, 3) + 500)) + nan;
        blank_frames(:, :, 1:size(stack, 3)) = stack;
        stack = blank_frames;
        clear blank_frames;
    else
    end
    
    try
        frame = imread(fpath, frame_n);
        stack(:, :, frame_n) = frame;
    catch
        %no more frames in tif stack
        disp([int2str(frame_n - 1) ' frames loaded.'])
        break
    end
    
end
stack(:, :, (frame_n):size(stack, 3) ) = [];


end