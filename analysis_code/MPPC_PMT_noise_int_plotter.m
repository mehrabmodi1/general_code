clear all
close all

paths = [{'C:\Data\Data\Raw_data\20190506\MPPC_pollen_grain_power_range\'}, ...
            {'C:\Data\Data\Raw_data\20190506\MPPC_pollen_grain_power_range_detectors_swapped\'}];

path = paths{1};        
cd(path);
MPPC_colour = [0.9, 0.4, 0.6];
PMT_colour = [0.5, 0.6, 0.8];

if exist([path, 'data_mat.mat']) ~= 2
    tif_list = dir('*.tif');
    ROI_mat = load([path, 'ROI_mat.mat']);
    ROI_mat = ROI_mat.ROI_mat;
    curr_pix = find(ROI_mat == 1);
    n_points = size(tif_list, 1);
    l_power_vec = [];

    for point_n = 1:n_points
        curr_name = tif_list(point_n).name;

        %parsing filename to get current laser power
        dashi = findstr(curr_name, '_');
        l_power_vec = [l_power_vec, str2num(curr_name(1:(dashi - 1)))];

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

        curr_pix_PMT = [];
        curr_pix_MPPC = [];
        for frame_n = 1:size(stack_PMT, 3)
            %computing mean and sd of intensity of lists of pixel vals
            curr_fr_PMT = stack_PMT(:, :, frame_n);
            curr_pix_PMT = [curr_pix_PMT; curr_fr_PMT(curr_pix)];

            curr_fr_MPPC = stack_MPPC(:, :, frame_n);
            curr_pix_MPPC = [curr_pix_MPPC; curr_fr_MPPC(curr_pix)];
        end
        
        mean_PMT_vec(point_n, 1) = mean(curr_pix_PMT);
        sd_PMT_vec(point_n, 1) = std(curr_pix_PMT);

        mean_MPPC_vec(point_n, 1) = mean(curr_pix_MPPC);
        sd_MPPC_vec(point_n, 1) = std(curr_pix_MPPC);


    end
    data_mat = [l_power_vec', mean_PMT_vec, sd_PMT_vec, mean_MPPC_vec, sd_MPPC_vec];
    save([path, 'data_mat.mat'], 'data_mat');

else
    data_mat = load([path, 'data_mat.mat']);
    data_mat = data_mat.data_mat;
    l_power_vec = data_mat(:, 1)';
    mean_PMT_vec = data_mat(:, 2);
    sd_PMT_vec = data_mat(:, 3);
    mean_MPPC_vec = data_mat(:, 4);
    sd_MPPC_vec = data_mat(:, 5);
end

%subtracting baseline obtained at lowest, ~ no illumination
mean_PMT_vec = mean_PMT_vec - repmat(mean_PMT_vec(1), size(mean_PMT_vec, 1), 1);
mean_MPPC_vec = mean_MPPC_vec - repmat(mean_MPPC_vec(1), size(mean_MPPC_vec, 1), 1);

norm_PMT_vec = mean_PMT_vec./sd_PMT_vec;
norm_MPPC_vec = mean_MPPC_vec./sd_MPPC_vec;

%computing some multiple of green light intensity to plot on x
Gintensity_vec = l_power_vec.^2;


plot(Gintensity_vec, norm_PMT_vec, '.', 'markerSize', 20, 'Color', PMT_colour)
hold on
plot(Gintensity_vec.*1.32, norm_MPPC_vec, '.', 'markerSize', 20, 'Color', MPPC_colour)      %multiplying x axis by 0.75 to correct for detection arm assymetry
ylabel('detector signal (mean/sd)')
xlabel('light intensity (AU)')
fig_wrapup(1, '')

figure(2)
loglog(Gintensity_vec, norm_PMT_vec, '.', 'markerSize', 20, 'Color', PMT_colour)
hold on
loglog(Gintensity_vec.*1.32, norm_MPPC_vec, '.', 'markerSize', 20, 'Color', MPPC_colour)          %multiplying x axis by 0.75 to correct for detection arm assymetry
ylabel('detector signal (mean/sd)')
xlabel('light intensity (AU)')
fig_wrapup(2, '')

figure(3)
plot(norm_PMT_vec, norm_MPPC_vec, '.', 'markerSize', 20, 'Color', 'k')
xlabel('PMT signal (mean/sd)')
ylabel('MPPC signal (mean/sd)')
fig_wrapup(3, [])


figure(4)
loglog(norm_PMT_vec, norm_MPPC_vec, '.', 'markerSize', 20, 'Color', 'k')
xlabel('PMT signal (mean/sd)')
ylabel('MPPC signal (mean/sd)')
fig_wrapup(4, [])


%comparing PMT mean F across detector positions
data_mat1 = load([paths{1}, 'data_mat.mat']);
data_mat1 = data_mat1.data_mat;
laser_vec1 = data_mat1(:, 1);
norm_PMT_vec1 = data_mat1(:, 2);

data_mat2 = load([paths{2}, 'data_mat.mat']);
data_mat2 = data_mat2.data_mat;
laser_vec2 = data_mat2(:, 1);
norm_PMT_vec2 = data_mat2(:, 2);

%plotting PMT norm means measured at the same laser value in each
%orientation
[common_l_vals, common_l_vals1i, common_l_vals2i] = intersect(laser_vec1, laser_vec2);

figure(5)
plot(norm_PMT_vec1(common_l_vals1i), norm_PMT_vec2(common_l_vals2i), '.', 'markerSize', 20, 'Color', 'k')
xlabel('mean F on PMT')
ylabel('mean F on PMT, position swapped')
axis([0, 2500, 0, 2500])
fig_wrapup(5, [])
pbaspect([1 1 1])
PMT_ratio = norm_PMT_vec1(common_l_vals1i)./norm_PMT_vec2(common_l_vals2i);
PMT_ratio = mean(PMT_ratio(7:end));