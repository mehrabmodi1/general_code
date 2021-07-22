clear all
close all

tif_path = 'C:\Data\Data\Analysed_data\Suite2P_results\20210413\fly1_d5HT1bopGC6f_simple_handover\1\file_00001.tif';

script_name = mfilename;

stim_frs = [101, 152];

stack_obj = ScanImageTiffReader(tif_path);
stack = stack_obj.data();
stack = permute(stack,[2 1 3]);

PMT_stack = stack(:, :, 1:2:end);
PMT_stack = align_image_rows(PMT_stack, []);

%manually measured background values
bk_val_PMT = 110;
%subtracting background
PMT_stack = PMT_stack - bk_val_PMT;

%mean frames
ave_im_PMT = mean(PMT_stack, 3);

%baseline frames
base_fr_PMT = mean(PMT_stack(:, :, 1:(stim_frs(1) - 10)), 3, 'omitnan');
base_fr_PMT_bar = add_scale_bar(base_fr_PMT, 0.18, 10, 1);
base_fr_PMT(base_fr_PMT < 70) = 0;                              %thresholding away bacground to hide it in dF/F image

figure(1)
PMT_im = real(base_fr_PMT_bar.^0.5);
imagesc(PMT_im, [0, max(max(PMT_im)).*0.8])
colormap('gray')
pbaspect([1 1 1])
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])


%resp period frames
resp_fr_PMT = mean(PMT_stack(:, :, stim_frs(1):stim_frs(2)), 3, 'omitnan');

%dF/F frames
dFF_fr_PMT = (resp_fr_PMT - base_fr_PMT)./base_fr_PMT;
dFF_fr_PMT(isinf(dFF_fr_PMT)) = 0;

figure(4)
imagesc(dFF_fr_PMT, [0, 2.5])
colormap('gray')
pbaspect([1 1 1])
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])
%colorbar
pbaspect([1 1 1])

highlight_resp_pix(3, stack, stim_frs, 1, 0.099, 0);

