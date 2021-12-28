clear all
close all

path = 'E:\Data\Raw_Data_Current\Resonant\20211210\fly1_zstck_NO68_2cvslp_00001\';
%stack_name = 'fly1_zstck_NO68_1cvslp_hipower_post_00001.tif';
stack_name = 'iamge.tif';

curr_path = [path, stack_name];
% 
% im_obj = ScanImageTiffReader(curr_path);
% [frame_time, zoom, n_chans] = SI_tif_info(im_obj);
% stack = im_obj.data();
% stack = permute(stack,[2 1 3]);
% 
% %splitting R and G channels if needed
% if n_chans == 2
%     stack1 = stack(:, :, 1:2:end);
%     stack2 = stack(:, :, 2:2:end);
% else
% end

im = imread(curr_path);
[im_a, lag] = align_image_rows(im, 12);
im_a = mat2gray(im_a);
imwrite(im2uint16(im_a(:, 13:end)), [path, 'iamge_aligned.tif'])
keyboard