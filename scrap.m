close all
clear all

%stack_obj = ScanImageTiffReader('C:\Users\Mehrab\Dropbox (HHMI)\Kaspar-Mehrab\Fseries_images_code\MPPC_pollen_grain_power_range\22_percent_00001.tif');
stack_obj = ScanImageTiffReader('C:\Users\Mehrab\Dropbox (HHMI)\Kaspar-Mehrab\Fseries_images_code\MPPC_pollen_grain_power_range\04_percent_00001.tif');
stack = stack_obj.data();

PMT_stack = double(stack(:, :, 1:2:end));
MPPC_stack = double(stack(:, :, 2:2:end));


[PMT_hist, PMT_bins] = hist(reshape(PMT_stack, 1, []), 100);
[MPPC_hist, MPPC_bins] = hist(reshape(MPPC_stack, 1, []), 100);


a = double(PMT_stack < 0);
playStack(a, 0.03, 1)