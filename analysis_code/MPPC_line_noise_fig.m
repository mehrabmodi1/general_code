clear all
close all

path = 'C:\Data\Data\Raw_data\20190506\MPPC_pollen_grain_power_range\00_percent_laser_shuttered_00001.tif';

stack_obj = ScanImageTiffReader(path);
[frame_time, zoom, n_chans, PMT_offsets] = SI_tif_info(stack_obj);

MPPC_color = [0.9, 0.4, 0.6];
PMT_color = [0.5, 0.6, 0.8];

stk = stack_obj.data();
stk = permute(stk,[2 1 3]);
stk = double(stk);

PMT_stack = stk(:, :, 1:2:end);
MPPC_stack = stk(:, :, 2:2:end);

%normalising stacks
PMT_min = min(min(min(PMT_stack)));
PMT_max = max(max(max(PMT_stack)));
PMT_cutoff = 400;                       %determined manually
PMT_cutoff = (PMT_cutoff - PMT_min)./(PMT_max - PMT_min);
PMT_stack = (PMT_stack - PMT_min)./(PMT_max - PMT_min);

MPPC_min = min(min(min(MPPC_stack)));
MPPC_max = max(max(max(MPPC_stack)));
MPPC_cutoff = -560;                     %determined manually
MPPC_cutoff = (MPPC_cutoff - MPPC_min)./(MPPC_max - MPPC_min);
MPPC_stack = (MPPC_stack - MPPC_min)./(MPPC_max - MPPC_min);

figure(1)
imagesc(PMT_stack(:, :, 1), [0, PMT_cutoff])
colormap('gray')
colorbar;
figure(2)
PMT_cutoff_fr = im2bw(PMT_stack(:, :, 1), PMT_cutoff);
imagesc(PMT_cutoff_fr)
colormap('gray')

figure(3)
imagesc(MPPC_stack(:, :, 1), [0, MPPC_cutoff])
colormap('gray')
colorbar;
figure(4)
PMT_cutoff_fr = im2bw(PMT_stack(:, :, 1), PMT_cutoff);
imagesc(PMT_cutoff_fr)
colormap('gray')
colorbar;

PMT_pix_vals = reshape(PMT_stack, [], 1);
MPPC_pix_vals = reshape(MPPC_stack, [], 1);


bin_vals = [0:0.005:1];
PMT_hist = hist(PMT_pix_vals, bin_vals);
PMT_hist = (PMT_hist./sum(PMT_hist)).*100;
MPPC_hist = hist(MPPC_pix_vals, bin_vals);
MPPC_hist = (MPPC_hist./sum(MPPC_hist)).*100;

figure(5)
plot(bin_vals, PMT_hist, 'lineWidth', 2, 'Color', PMT_color)
hold on
plot([PMT_cutoff, PMT_cutoff], [0, 100], '--r', 'lineWidth', 2);
axis([0, 1, 0, 0.05])
fig_wrapup(5, [])
figure(6)
plot(bin_vals, MPPC_hist, 'lineWidth', 2, 'Color', MPPC_color)
hold on
plot([MPPC_cutoff, MPPC_cutoff], [0, 100], '--r', 'lineWidth', 2);
axis([0, 1, 0, 0.05])
fig_wrapup(6, [])


