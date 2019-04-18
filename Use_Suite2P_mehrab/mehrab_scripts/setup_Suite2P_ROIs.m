function [ROI_mat] = setup_Suite2P_ROIs(data_mat)

n_rois_all = size(data_mat.stat, 2);
list_rois = zeros(n_rois_all, 1);
for roi_n = 1:n_rois_all
    list_rois(roi_n, 1) = data_mat.stat(roi_n).iscell; 
end

list_cell_rois = find(list_rois == 1); 
n_rois = length(list_cell_rois);

size_im = [data_mat.cl.Ly, data_mat.cl.Lx];
ROI_mat = zeros(size_im(1), size_im(2), n_rois);        %each ROI separated along third dimension.
for roi_n = 1:n_rois
    roi_n_orig = list_cell_rois(roi_n);
    roi_pix = data_mat.stat(roi_n_orig).ipix;           %pixels for current roi
    ROI_im = zeros(size_im(1), size_im(2));
    ROI_im(roi_pix) = 1;
    ROI_mat(:, :, roi_n) = ROI_im;                      %Suite2P generates an ROI_mat that is smaller because of the alingment shifts
end

%making ROI mat larger to match the size of the registered frames for proper data re-extraction.
mean_im = data_mat.ops.mimg;        %original size of frames
base_mat = zeros(size(mean_im, 1), size(mean_im, 2), size(ROI_mat, 3));
base_mat(1:size(ROI_mat, 1), 1:size(ROI_mat, 2), 1:size(ROI_mat, 3)) = ROI_mat;
ROI_mat_orig = ROI_mat;
ROI_mat = base_mat;
clear base_mat

%aligning ROI matrix to underlying mean image
ROI_mat_squashed = max(ROI_mat, [], 3);
sub_im = data_mat.mimg(:, :, 2);
sub_im = sub_im - mean(mean(sub_im));
mean_im = mean_im - mean(mean(mean_im));
corr_mat = xcorr2(mean_im, ROI_mat_squashed);
[del, lags] = max(corr_mat(:));
[xlag, ylag] = ind2sub(size(corr_mat), lags);
xlag = xlag - size(mean_im, 1);
ylag = ylag - size(mean_im, 2);

ROI_mat1 = circshift(ROI_mat, xlag, 1);
ROI_mat1 = circshift(ROI_mat1, ylag, 2);
ROI_mat_squashed = max(ROI_mat1, [], 3);

% disp_im = mean_im + ROI_mat_squashed.*max(max(mean_im));
% subplot(2, 1, 1)
% imagesc(disp_im)

% subplot(2, 1, 2)
% disp_im2 = data_mat.mimg(:, :, 2) +  max(ROI_mat_orig, [], 3).*max(max(mean_im));
% imagesc(disp_im2)
% title('are ROIs correctly aligned to mean image?')

end
