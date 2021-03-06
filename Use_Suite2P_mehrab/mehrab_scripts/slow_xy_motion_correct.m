function [reg_stack, saved_lags, bk_ROI] = slow_xy_motion_correct(curr_stack, ref_im, lag_mat, detailed_lag_mat, registration_type)
%syntax: reg_stack = slow_xy_motion_correct(curr_stack, ref_im, lag_mat, registration_type)
%If registration_type = 1, this function assumes any x-y motion is slow and negligible within a
%single trial (tif stack). It then averages all frames in curr_stack and
%aligns the result to ref_im. It uses the same lags to correct the
%individual frames in curr_stack to give reg_stack.
%If registration_type = 2, this function first does the whole-stack
%correction and then computes lags for aach frame to do a per-frame motion
%correction.
%Mehrab Modi 20180118

%test variables easy
% ref_im = zeros(100, 100);
% ref_im(40:50, 40:50) = 5;
% ref_im = ref_im + 1;
% curr_stack = zeros(100, 100);
% curr_stack(30:40, 30:40) = 5;
% curr_stack = repmat(curr_stack, 1, 1, 30);
% curr_stack = curr_stack + 1;

% %test variables real
% stack1_path = 'C:\Data\Data\Raw_data\20180111\fly1_axons_train_stim\fly1_od_trains_00013.tif';
% stack2_path = 'C:\Data\Data\Raw_data\20180111\fly1_axons_train_stim\fly1_od_trains_00089.tif';
% stack1 = ScanImageTiffReader(stack1_path).data();
% curr_stack = ScanImageTiffReader(stack2_path).data();
% ref_im = mean(stack1, 3, 'omitnan');

%testing sign convention of lag mat
%lag_mat(:, 1:2) = lag_mat(:, 1:2) * -1;

if isempty(lag_mat) == 1
    curr_im = mean(curr_stack, 3, 'omitnan');
    c = xcorr2_fft(curr_im, ref_im); 
    [maxr, maxcol] = find(c == max(max(c)));

    %computing lags
    col_lag = size(curr_stack, 1) - maxr;
    row_lag = size(curr_stack, 2) - maxcol;

%case where lags for this trial have been determined by a manually selected
%landmark
elseif isempty(lag_mat) == 0
    col_lag = round(lag_mat(1, 2));
    row_lag = round(lag_mat(1, 1));
end
 
%generating a warning if lags more than 20% of size of frame
if mean([col_lag, row_lag]) > (size(ref_im, 1)./5)
    disp('WARNING: X-Y movement of more than 20% frame size detected.')
else
end

%shifting stack
reg_stack = translate_stack(curr_stack, [row_lag; col_lag], nan);

%doing per-frame registration, if specified for by user
if registration_type == 2
   reg_stack_orig = reg_stack;
   reg_stack = zeros(size(reg_stack_orig, 1), size(reg_stack_orig, 2), size(reg_stack_orig, 3)) + nan;
   ref_im = mean(reg_stack_orig, 3, 'omitnan');
   for frame_n = 1:size(reg_stack, 3)
       if isempty(detailed_lag_mat) == 1

           curr_im = reg_stack_orig(:, :, frame_n);
           c = xcorr2_fft(curr_im, ref_im); 
           [maxr, maxcol] = find(c == max(max(c)));

           %computing lags
           col_lag = size(curr_stack, 1) - maxr;
           row_lag = size(curr_stack, 2) - maxcol;

           if sign(col_lag) == 1
               col_lag = min([col_lag, 20]);
           elseif sign(col_lag) == -1
               col_lag = max([col_lag, -20]);
           else
           end

           if sign(row_lag) == 1
               row_lag = min([row_lag, 20]);
           elseif sign(row_lag) == -1
               row_lag = max([row_lag, -20]);
           else
           end
       else
           keyboard
       end
       
       saved_lags(frame_n, :) = [row_lag, col_lag];
       curr_im = translate_stack(curr_im, [row_lag; col_lag], nan);     %note: for int16, nan is not defined and set to 0 by matlab.
       reg_stack(:, :, frame_n) = curr_im;

   end
   
   
    
    
elseif registration_type == 3
    for frame_n = 1:size(curr_stack, 3)
        curr_im = curr_stack(:, :, frame_n);
        row_lag = detailed_lag_mat(frame_n, 1);
        col_lag = detailed_lag_mat(frame_n, 2);
        try
            curr_im = translate_stack(curr_im, [row_lag; col_lag], nan);
        catch
            keyboard
        end
        reg_stack(:, :, frame_n) = curr_im;
        saved_lags(frame_n, :) = [row_lag, col_lag];
       
    end
else
    saved_lags = [];
end






%testing plots
% ave_corrected = mean(reg_stack, 3, 'omitnan');
% ave_uncorrected = mean(curr_stack, 3, 'omitnan');
% figure(1)
% subplot(2, 2, 1)
% imagesc(ref_im)
% subplot(2, 2, 2)
% imagesc(ave_uncorrected)
% subplot(2, 2, 3)
% imagesc(ave_corrected)
% subplot(2, 2, 4)
% imagesc(ref_im - ave_corrected)
% 
% figure(2)
% plot(saved_lags)
% 
% del = input('press enter');
% close figure 1
% close figure 2
% try
%     close figure 3
% catch
% end
