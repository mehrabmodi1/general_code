clear all
close all

path = 'C:\Data\Data\Raw_data\20190327\MPPC_PMT_Fplatic\';
cd(path);
tif_list = dir('*.tif');
ROI = load([path, 'ROI.mat']);
ROI_mat = ROI.ROI;
curr_pix = find(ROI_mat == 1);
n_points = size(tif_list, 1);

MPPC_colour = [0.9, 0.4, 0.6];
PMT_colour = [0.5, 0.6, 0.8];

for point_n = 1:n_points
    curr_name = tif_list(point_n).name;
    
    %reading in stack object
    im_obj = ScanImageTiffReader([path, curr_name]);
    %obtaining image stack
    stack = im_obj.data();
    %stack_orig = stack;
    stack = permute(stack,[2 1 3]);
    stack = double(stack);
    
    %separating channels
    stack_PMT = stack(:, :, 1:2:size(stack, 3));
    stack_MPPC = stack(:, :, 2:2:size(stack, 3));
    
    %computing mean intensity and sd of intensity
    %lists of pixel vals
    curr_fr_PMT = stack_PMT(:, :, 50);
    curr_pix_PMT = curr_fr_PMT(curr_pix);
    
    curr_fr_MPPC = stack_MPPC(:, :, 50);
    curr_pix_MPPC = curr_fr_MPPC(curr_pix);
    
    mean_PMT_vec(point_n, 1) = mean(curr_pix_PMT);
    sd_PMT_vec(point_n, 1) = std(curr_pix_PMT);
    
    mean_MPPC_vec(point_n, 1) = mean(curr_pix_MPPC);
    sd_MPPC_vec(point_n, 1) = std(curr_pix_MPPC);
    
   
end

%subtracting baseline obtained at lowest, ~ no illumination
mean_PMT_vec = mean_PMT_vec - repmat(mean_PMT_vec(1), size(mean_PMT_vec, 1), 1);
mean_MPPC_vec = mean_MPPC_vec - repmat(mean_MPPC_vec(1), size(mean_MPPC_vec, 1), 1);

norm_PMT_vec = mean_PMT_vec./sd_PMT_vec;
norm_MPPC_vec = mean_MPPC_vec./sd_MPPC_vec;

%computing some multiple of green light intensity to plot on x
l_power_vec = [0, 0.5, 1, 2, 4, 6, 8, 10];
Gintensity_vec = l_power_vec.^2;


plot(Gintensity_vec, norm_PMT_vec, 'Color', PMT_colour, 'LineWidth', 2)
hold on
plot(Gintensity_vec, norm_MPPC_vec, 'Color', MPPC_colour, 'LineWidth', 2)
ylabel('detector signal (mean/sd)')
xlabel('light intensity (AU)')
fig_wrapup(1, '')

figure(2)
loglog(Gintensity_vec, norm_PMT_vec, 'Color', PMT_colour, 'LineWidth', 2)
hold on
loglog(Gintensity_vec, norm_MPPC_vec, 'Color', MPPC_colour, 'LineWidth', 2)
ylabel('detector signal (mean/sd)')
xlabel('light intensity (AU)')
fig_wrapup(2, '')

keyboard

PMT_pkresp_val = (121.3 + 113.5);
MPPC_pkresp_val = (30.02 + 290.8);
