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
    var_PMT_vec(point_n, 1) = std(curr_pix_PMT).^2;
    
    mean_MPPC_vec(point_n, 1) = mean(curr_pix_MPPC);
    var_MPPC_vec(point_n, 1) = std(curr_pix_MPPC).^2;
    
    
end

%plotting
mean_PMT_vec = mean_PMT_vec;
mean_MPPC_vec = mean_MPPC_vec - min(mean_MPPC_vec);
val_vecs_PMT = [mean_PMT_vec, var_PMT_vec];
val_vecs_PMT = sortrows(val_vecs_PMT);
val_vecs_MPPC = [mean_MPPC_vec, var_MPPC_vec];
val_vecs_MPPC = sortrows(val_vecs_MPPC);
val_vecs_PMT = [val_vecs_PMT, val_vecs_PMT(:, 2).^0.5];
val_vecs_MPPC = [val_vecs_MPPC, val_vecs_MPPC(:, 2).^0.5];

PMT_pkresp_val = (121.3 + 113.5);
MPPC_pkresp_val = (30.02 + 290.8);

figure(1)
plot(val_vecs_PMT(:, 1), val_vecs_PMT(:, 2), 'o', 'Color', PMT_colour, 'MarkerSize', 6, 'MarkerFaceColor', PMT_colour);
hold on
plot(val_vecs_MPPC(:, 1), val_vecs_MPPC(:, 2), 'o', 'Color', MPPC_colour, 'MarkerSize', 6, 'MarkerFaceColor', MPPC_colour);
plot(val_vecs_PMT(:, 1), val_vecs_PMT(:, 2), 'Color', PMT_colour);
plot(val_vecs_MPPC(:, 1), val_vecs_MPPC(:, 2), 'Color', MPPC_colour);
xlabel('mean F (AU)');
ylabel('variance F (AU^2)');
script_name = mfilename;
fig_wrapup(1, script_name)

figure(2)
shadedErrorBar([], val_vecs_PMT(:, 1), val_vecs_PMT(:, 3),{'-o', 'Color', PMT_colour, 'MarkerSize', 6, 'MarkerFaceColor', PMT_colour});
hold on
shadedErrorBar([], val_vecs_MPPC(:, 1), val_vecs_MPPC(:, 3),{'-o', 'Color', MPPC_colour, 'MarkerSize', 6, 'MarkerFaceColor', MPPC_colour});
xlabel('illumination intensity');
ylabel('mean F +/- SD (AU)')
fig_wrapup_nonums(2, script_name, [1]);
axis([0, 8, -2000, 10000]);