clear all
close all

path = 'C:\Data\Data\Analysed_data\Manual_ROI_results\20190312\MPPC_PMT_13F02LexA_opGC6f_short_dur - Copy\';
%reading in stack object
im_obj = ScanImageTiffReader([path, 'odor_trs_00001.tif']);
%obtaining image stack
stack = im_obj.data();
%stack_orig = stack;
stack = permute(stack,[2 1 3]);
stack = double(stack);

%separating channels
stack_PMT = stack(:, :, 1:2:size(stack, 3));
stack_MPPC = stack(:, :, 2:2:size(stack, 3));

bk_ROI = load([path, 'bk_ROI.mat']);
bk_ROI = bk_ROI.bk_ROI;
curr_pix = find(bk_ROI == 1);

curr_fr_PMT = stack_PMT(:, :, 792);
bk_val = mean(curr_fr_PMT(curr_pix))

curr_fr_MPPC = stack_MPPC(:, :, 792);
bk_val = mean(curr_fr_MPPC(curr_pix))
