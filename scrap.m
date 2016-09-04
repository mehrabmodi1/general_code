clear all
close all

path = 'D:\Data\CSHL\Analysed_data\Marius\Results\mouse_name\date\1\F_mouse_name_date_plane1_Nk650';

all_data = load(path);
ROI_data = all_data.stat;

n_rois = size(ROI_data, 2);
ROI_mat = zeros(256, 255);
for roi_n = 1:n_rois
    roi_pixi = ROI_data(roi_n).ipix; 
    ROI_mat(roi_pixi) = roi_n;    
end

figure(1)
imagesc(ROI_mat)