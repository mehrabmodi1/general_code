clear all
close all

dir = 'C:\Data\Data\Analysed_data\Analysis_results\subtract_linenoise\';
stack_obj = ScanImageTiffReader([dir, 'odor_trs_00013.tif']); 
stack = stack_obj.data();
stack = permute(stack, [2, 1, 3]);

n_frames = size(stack, 3);
eg_lin = medfilt(stack());

%%PICK UP THREAD HERE
%need to figure out how to subtract away line noise from image background.
%maybe 1D median filter the whole frame or try out 2D median filt.

