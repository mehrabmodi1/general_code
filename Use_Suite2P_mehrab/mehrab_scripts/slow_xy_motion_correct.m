function reg_stack = slow_xy_motion_correct(curr_stack, ref_im)
%syntax: reg_stack = slow_xy_motion_correct(curr_stack, ref_im)
%This function assumes any x-y motion is slow and negligible within a
%single trial (tif stack). It then averages all frames in curr_stack and
%aligns the result to ref_im. It uses the same lags to correct the
%individual frames in curr_stack to give reg_stack.
%Mehrab Modi 20180118

% %test variables
% stack1_path = 'C:\Data\Data\Raw_data\20180111\fly1_axons_train_stim\fly1_od_trains_00013.tif';
% stack2_path = 'C:\Data\Data\Raw_data\20180111\fly1_axons_train_stim\fly1_od_trains_00089.tif';
% stack1 = ScanImageTiffReader(stack1_path).data();
% curr_stack = ScanImageTiffReader(stack2_path).data();
% curr_stack_orig = curr_stack;
% ref_im = mean(stack1, 3, 'omitnan');

curr_im = mean(curr_stack, 3, 'omitnan');
c = xcorr2_fft(curr_im, ref_im); 

[maxr, maxcol] = find(c == max(max(c)));

%computing lags
col_lag = size(curr_stack, 1) - maxr;
row_lag = size(curr_stack, 2) - maxcol;

%shifting stack
reg_stack = circshift(curr_stack, row_lag, 1);
reg_stack = circshift(reg_stack, col_lag, 2);

%replacing circularly shifted pixels with nans
if sign(col_lag) == 1
    reg_stack(:, 1:col_lag, :) = nan;
elseif sign(col_lag) == -1
    reg_stack(:, (end + col_lag):end, :) = nan;
else
end

if sign(row_lag) == 1
    reg_stack(1:row_lag, :, :) = nan;
elseif sign(row_lag) == -1
    reg_stack((end + row_lag):end, :, :) = nan;
else
end


% %testing plots
% ave_corrected = mean(reg_stack, 3, 'omitnan');
% figure(1)
% subplot(3, 1, 1)
% imagesc(ref_im)
% subplot(3, 1, 2)
% imagesc(ave_corrected)
% subplot(3, 1, 3)
% imagesc(ref_im - ave_corrected)

