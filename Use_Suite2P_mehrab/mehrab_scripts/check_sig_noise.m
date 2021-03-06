function [noise_mat] = check_sig_noise(frame)
%Syntax:[noise_mat] = check_sig_noise(frame)
%This function uses columns in frame to compute estimates of slowly changing 
%noise (such as line noise) in the PMT signal. Any slow noise should be unchaging
%across pixels in a row. 



%picking vertical column ROIs at lateral edges to estimate background noise
%contribution in each pixel row without contribution from biological sample.
row_length = size(frame, 2);
col_width = round(row_length./10);

left_bk_vals = mean(frame(:, 1:col_width), 2);
left_bk_vals = movmean(left_bk_vals, 10);
right_bk_vals = mean(frame(:, (row_length - col_width):row_length), 2);
right_bk_vals = movmean(right_bk_vals, 10);

c = corrcoef(left_bk_vals, right_bk_vals);

if c(1, 2) > 0.8
    noise_vec = mean([left_bk_vals, right_bk_vals], 2);
    noise_mat = repmat(noise_vec, 1, row_length);
    noise_mat = double(noise_mat);
    
else
    noise_mat = zeros(size(frame, 1), size(frame, 2));
    noise_mat = double(noise_mat);
    
end

%function testing lines
% imagesc(frame)
% figure(2)
% imagesc(left_bk_vals)
% title('left edge mean vals')
% figure(3)
% imagesc(right_bk_vals)
% title('right edge mean vals')
% figure(4)
% corr_frame = frame - noise_mat;
% imagesc(corr_frame)

