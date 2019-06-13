clear all
close all

tif_path = 'C:\Data\Data\Raw_data\20190506\fly2_MPPC_OK107_opGC6f_highpower\odor_trs7.4percent_00003.tif';

script_name = mfilename;

stack_obj = ScanImageTiffReader(tif_path);
stack = stack_obj.data();
stack = permute(stack,[2 1 3]);

PMT_stack = stack(:, :, 1:2:end);
MPPC_stack = stack(:, :, 2:2:end);
PMT_stack = align_image_rows(PMT_stack, []);
MPPC_stack = align_image_rows(MPPC_stack, 1);

%manually measured background values
bk_val_MPPC = -650;
bk_val_PMT = 110;
%subtracting background
PMT_stack = PMT_stack - bk_val_PMT;
MPPC_stack = MPPC_stack - bk_val_MPPC;

%mean frames
ave_im_PMT = mean(PMT_stack, 3);
ave_im_MPPC = mean(MPPC_stack, 3);


%baseline frames
base_fr_PMT = mean(PMT_stack(:, :, 1:310), 3, 'omitnan');
base_fr_PMT_bar = add_scale_bar(base_fr_PMT, 0.18, 10, 1);
base_fr_PMT(base_fr_PMT < 70) = 0;                              %thresholding away bacground to hide it in dF/F image
base_fr_MPPC = mean(MPPC_stack(:, :, 1:310), 3, 'omitnan');
base_fr_MPPC_bar = add_scale_bar(base_fr_MPPC, 0.18, 10, 1);
base_fr_MPPC(base_fr_MPPC < 7) = 0;                             %thresholding away bacground to hide it in dF/F image

figure(1)
PMT_im = real(base_fr_PMT_bar.^0.5);
imagesc(PMT_im, [0, max(max(PMT_im)).*0.8])
colormap('gray')
pbaspect([1 1 1])
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])


figure(2)
MPPC_im = real(base_fr_MPPC_bar.^0.5);
imagesc(MPPC_im, [0, max(max(MPPC_im)).*0.8])
colormap('gray')
pbaspect([1 1 1])
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])




%resp period frames
resp_fr_PMT = mean(PMT_stack(:, :, 341:475), 3, 'omitnan');
resp_fr_MPPC = mean(MPPC_stack(:, :, 341:475), 3, 'omitnan');

%dF/F frames
dFF_fr_PMT = (resp_fr_PMT - base_fr_PMT)./base_fr_PMT;
dFF_fr_PMT(isinf(dFF_fr_PMT)) = 0;
dFF_fr_MPPC = (resp_fr_MPPC - base_fr_MPPC)./base_fr_MPPC;
dFF_fr_MPPC(isinf(dFF_fr_MPPC)) = 0;

figure(3)
imagesc(dFF_fr_PMT, [0, 2.5])
colormap('gray')
pbaspect([1 1 1])
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])
%colorbar
pbaspect([1 1 1])

figure(4)
imagesc(dFF_fr_MPPC, [0, 2.5])
colormap('gray')
pbaspect([1 1 1])
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])
%colorbar
pbaspect([1 1 1])
