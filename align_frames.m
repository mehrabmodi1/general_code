function [stack_r, lags] = align_frames(varargin)
%[stack_r, lags] = align_frames(stack, target, thresh);
%Function to align the frames in the 3d matrix stack to the image in
%target. Alignment based on rigid translation (no rotation) according to
%lags in a 2D cross-correlation. Thresh is a boolean variable, default 1.
%If 1, the function thresholds target and each frame before calculating 2D
%xcorr. This runs quicker and should be less affected by activity related
%intensity changes.


% direc = 'C:\Data\CSHL\test direc\';
% stack = zeros(256, 256, 100);
% for frame_n = 1:100
%     stack(:, :, frame_n) = (imread([direc 'test_direc.tif'], frame_n) );
%     
% end
% target = stack(:, :, 1);

if nargin == 3
     stack = varargin{1};
     target = varargin{2};
     thresh = varargin{3};
elseif nargin == 2 
    stack = varargin{1};
    n_frames = size(stack, 3);
    target = varargin{2};
    thresh = 1;
elseif nargin == 1
    stack = varargin{1};
    n_frames = size(stack, 3);
    n_ave = floor(n_frames./ 10);
    target = mean(stack(:, :, 1:n_ave), 3);
    thresh = 1;
else
end
    
if thresh == 1
    %thresholding and binarising target image
    target = im2bw(target, graythresh(target));
    %getting rid of small objects
    target = single(bwareaopen(target, 150));
    
else
end
 
stack_r = zeros(size(stack, 1), size(stack, 2), size(stack, 3)) + nan;
lags = zeros(n_frames, 2);
for frame_n = 1:n_frames
    frame = stack(:, :, frame_n);
    
    %identifying lags using 2d cross correlation
    
    if thresh == 1
        %thresholding and binarising target image
        frame_BW = im2bw(frame, graythresh(frame));
        %getting rid of small objects
        frame_BW = single(bwareaopen(frame_BW, 150));
    else
    end
    
    lag_mat = xcorr2(target, frame_BW);                %generating the 2D cross-correlation matrix
    max_corr = max(max(lag_mat));                   %identifying the highest correlation at any lag values
    [row_lag, col_lag] = find(lag_mat == max_corr);     %identifying x and y lags of max corr
    row_lag = row_lag - size(target, 1);                            %xcorr mat has double the width and height of target/stack
    col_lag = col_lag - size(target, 2);
    row_lag = min(row_lag);
    col_lag = min(col_lag);
    
    clear lag_mat
    clear max_corr
    
    %offsetting current frame and padding missing area with nans.
    frame_r = zeros(size(frame, 1), size(frame, 2)) + nan;
    n_rows = size(frame, 1);
    n_cols = size(frame, 2);
    diff_r = abs(row_lag);
    diff_c = abs(col_lag);
    
    sign_r = sign(row_lag);
    if sign_r == 0
        sign_r = 1;
    else
    end
    
    sign_c = sign(col_lag);
    if sign_c == 0
        sign_c = 1;
    else
    end
    
    
    if sign_r > 0 && sign_c > 0                                                   %Case1: frame needs to move downward and rightward
        frame_seg = frame( (1):(n_rows - diff_r), (1):(n_cols - diff_c) );
        frame_r( (diff_r + 1):(n_rows),(diff_c + 1):(n_cols) ) = frame_seg;     %inserting segment of frame into registered frame
        
    elseif sign_r < 0 && sign_c > 0                                               %Case2: frame needs to move upward and rightward
        frame_seg = frame( (diff_r + 1):(n_rows), (1):(n_cols - diff_c) );      
        frame_r( (1):(n_rows - diff_r),(diff_c + 1):(n_cols) ) = frame_seg;     %inserting segment of frame into registered frame
        
    elseif sign_r > 0 && sign_c < 0                                               %Case3: frame needs to move downward and leftward
        frame_seg = frame( (1):(n_rows - diff_r), (diff_c + 1):(n_cols) );
        frame_r( (diff_r + 1):(n_rows),(1):(n_cols - diff_c) ) = frame_seg;     %inserting segment of frame into registered frame
        
    elseif sign_r < 0 && sign_c < 0                                               %Case4: frame needs to move upward and leftward
        frame_seg = frame( (diff_r + 1):(n_rows), (diff_c + 1):(n_cols) );
        frame_r( (1):(n_rows - diff_r),(1):(n_cols - diff_c) ) = frame_seg;     %inserting segment of frame into registered frame
    
    end
    
    
    stack_r(:, :, frame_n) = frame_r;
    lags = [row_lag, col_lag];
end

end  


