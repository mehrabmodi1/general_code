clear all
close all

%base_dir = 'C:\Data\Data\Adithya_Yarena_data\5sLED\center_on\';
base_dir = 'C:\Data\Data\Adithya_Yarena_data\5sLED\periph_on\';

dir_contents = dir(base_dir);
dir_contents(1:2) = [];
traj_win = 5;   %time window in s sampled for plotting
frame_time = 2./11;
traj_win = round(traj_win./frame_time);

density_mat_e = zeros(1280, 1024);

%occupancy score matrix paramters
n_interp = 50;  %number of points to interpolate in between every pair of frames
blur_size = 2;  %SD of gaussian used to blur matrix

arena_ROI = load([base_dir, 'arena_ROI.mat']);
arena_ROI = arena_ROI.arena_ROI;

ROI_boundaries = [303, 682, 362, 746; 396, 501, 483, 501; 586, 677, 550, 756];

traj_samples_all = [];
saved_scores = [];

for dir_n = 1:size(dir_contents, 1)
    if isdir([base_dir, dir_contents(dir_n).name, '\']) == 0
        continue
    else
    end
    curr_dir = [base_dir, dir_contents(dir_n).name, '\'];
    
    data_struc = load([curr_dir, 'all_variables.mat']);
        
    coords = data_struc.xy;
    led_vec = data_struc.led_vec;
    
    %cleaning up tracking data
    del = find(coords(1, :) == -1 & coords(2, :) == -1);
    coords(:, del) = [];
    led_vec(del) = [];
    
    %getting stimulus defining ROIs
    masks = zeros(1024, 1280, 4);
    for maski = 1:4
        eval(['curr_mask = data_struc.binarymask', num2str(maski), ';']);
        del = find(arena_ROI' == 1 & curr_mask == 1);
        curr_mask = zeros(1024, 1280);
        curr_mask(del) = 1;
        eval(['masks(:, :, ', num2str(maski), ') = curr_mask;']);
               
    end
    
    
    %computing ROIs for occupancy quantifications
    se = strel('disk', 12);
    arena_roi_er = imerode(arena_ROI, se);     %eroding to exclude high counts ar edges of arena
        
    %original rois
    periph_roi = sum(masks(:, :, 2:4), 3)';
    center_roi = arena_ROI - periph_roi;
    
    %dilated rois
    se_dil = strel('disk', 10);
    periph_roi_dil = periph_roi;
    center_roi_dil = center_roi;
    for n_dilations = 1:5
        periph_roi_dil = imdilate(periph_roi_dil, se_dil);
        center_roi_dil = imdilate(center_roi_dil, se_dil);
    end
    
    quant_band = zeros(1024, 1280)';
    del = find(periph_roi_dil == 1 & center_roi_dil == 1);
    quant_band(del) = 1;
    
    %quantification rois
    periph_roi_quant = zeros(1280, 1024);
    del = find(quant_band == 1 & periph_roi == 1);
    periph_roi_quant(del) = 1;
    center_roi_quant = zeros(1280, 1024);
    del = find(quant_band == 1 & center_roi == 1);
    center_roi_quant(del) = 1;
    
    %excluding arena edges
    periph_roi_quant_er = zeros(1280, 1024);
    del = find(periph_roi_quant == 1 & arena_roi_er == 1);
    periph_roi_quant_er(del) = 1;
    center_roi_quant_er = zeros(1280, 1024);
    del = find(center_roi_quant == 1 & arena_roi_er == 1);
    center_roi_quant_er(del) = 1;
    
    
    %sampling traj segments around LED pulses
    d_led_vec = diff(led_vec);
    on_i = find(d_led_vec == 1);
    off_i = find(d_led_vec == -1);
    
    %accounting for interrupted LED pulse
    if length(off_i) == length(on_i) - 1
        on_i(end) = [];
    else
    end
    pulse_dur = median(off_i - on_i);       %in frames
    traj_samples = [];
    for pulse_n = 1:size(on_i, 2)
        if (on_i(pulse_n) - traj_win) < 0 || (off_i(pulse_n) + traj_win) > size(coords, 2)
            continue
        else
        end
        traj_samp = coords(:, (on_i(pulse_n) - traj_win):(off_i(pulse_n) + traj_win));
        traj_samples = pad_n_concatenate(traj_samples, coords(:, (on_i(pulse_n) - traj_win):(off_i(pulse_n) + traj_win)), 1, nan);
        
    end
    
    %integrating into location density matrix for entire trajectory
    for frame_n = 1:(size(coords, 2) - 1)
        curr_coords = coords(:, frame_n:(frame_n + 1) );
        curr_coords = round([linspace(curr_coords(1, 1), curr_coords(1, 2), n_interp); linspace(curr_coords(2, 1), curr_coords(2, 2), n_interp)]);     %interpolating between current pair of coords
    
        for interp_pt = 1:n_interp
            density_mat_e(curr_coords(1, interp_pt), curr_coords(2, interp_pt)) = density_mat_e(curr_coords(1, interp_pt), curr_coords(2, interp_pt)) + 1;        
        end
    end
    
    %quantifying density and scoring for current fly
    density_mat_ei = generate_occ_map(coords, n_interp).*arena_roi_er;
    density_map_off = generate_occ_map(traj_samples(:, (traj_win + pulse_dur):end), n_interp).*arena_roi_er;
    
    %scoring
    center_score_entire = mean(mean(density_mat_ei.*center_roi_quant_er));
    periph_score_entire = mean(mean(density_mat_ei.*periph_roi_quant_er));
    entire_score_er = periph_score_entire./center_score_entire;
    center_score_off = mean(mean(density_map_off.*center_roi_quant_er));
    periph_score_off = mean(mean(density_map_off.*periph_roi_quant_er));
    off_score_er = periph_score_off./center_score_off;
    saved_scores = [saved_scores; entire_score_er, off_score_er];
        
    traj_samples_all = pad_n_concatenate(traj_samples_all, traj_samples, 1, nan);
end

%play_trajs(traj_samples, masks, traj_win, pulse_dur, 0.075)

%plotting density map for entire trajectories
figure(1)
density_mat_e_er = density_mat_e.*arena_roi_er;
density_mat_e_bl = imgaussfilt(density_mat_e_er, blur_size);
set(gcf, 'Name', 'occupancy score plot, entire trajectories')
imagesc(density_mat_e_bl', [0, 10])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
colormap(gray)
hold on
plot(ROI_boundaries(1, [1, 3]), ROI_boundaries(1, [2, 4]), 'r', 'lineWidth', 2)
plot(ROI_boundaries(2, [1, 3]), ROI_boundaries(2, [2, 4]), 'r', 'lineWidth', 2)
plot(ROI_boundaries(3, [1, 3]), ROI_boundaries(3, [2, 4]), 'r', 'lineWidth', 2)
hold off
ax_vals = axis;
ax_vals(2) = ax_vals(2) - 400;
ax_vals(3) = ax_vals(3) + 130;
ax_vals(4) = ax_vals(4) - 85; 
axis(ax_vals);
fig_wrapup(1, [])

%getting density matrix for trajectories sampled on both sides of
%LED stim
density_map_s = generate_occ_map(traj_samples_all, n_interp).*arena_roi_er;
density_map_s_bl = imgaussfilt(density_map_s, blur_size);

figure(2)
set(gcf, 'Name', 'occupancy score plot, sampled trajectories')
imagesc(density_map_s_bl', [0, 1])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
colormap(gray)
hold on
plot(ROI_boundaries(1, [1, 3]), ROI_boundaries(1, [2, 4]), 'r', 'lineWidth', 2)
plot(ROI_boundaries(2, [1, 3]), ROI_boundaries(2, [2, 4]), 'r', 'lineWidth', 2)
plot(ROI_boundaries(3, [1, 3]), ROI_boundaries(3, [2, 4]), 'r', 'lineWidth', 2)
hold off
axis(ax_vals);
fig_wrapup(2, [])

%getting density matrix for trajectories sampled after LED stim ends
density_map_off = generate_occ_map(traj_samples_all(:, (traj_win + pulse_dur):end), n_interp).*arena_roi_er;
density_map_off_bl = imgaussfilt(density_map_off, blur_size);
figure(3)
set(gcf, 'Name', 'occupancy score plot, sampled after LED off')
imagesc(density_map_off_bl', [0, 1])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
colormap(gray)
hold on
plot(ROI_boundaries(1, [1, 3]), ROI_boundaries(1, [2, 4]), 'r', 'lineWidth', 2)
plot(ROI_boundaries(2, [1, 3]), ROI_boundaries(2, [2, 4]), 'r', 'lineWidth', 2)
plot(ROI_boundaries(3, [1, 3]), ROI_boundaries(3, [2, 4]), 'r', 'lineWidth', 2)
hold off
axis(ax_vals);
fig_wrapup(3, [])

%plotting trajectories sampled during LED-on period
figure(4)
plot(traj_samples_all(1:2:end, 1:traj_win), traj_samples_all(2:2:end, 1:traj_win), '.', 'Color', [0.65, 0.65, 0.65]);
hold on
plot(traj_samples_all(1:2:end, traj_win:pulse_dur), traj_samples_all(2:2:end, traj_win:pulse_dur), '.r');
plot(traj_samples_all(1:2:end, (traj_win + pulse_dur):end), traj_samples_all(2:2:end, (traj_win + pulse_dur):end), '.k');
plot(ROI_boundaries(1, [1, 3]), ROI_boundaries(1, [2, 4]), 'r', 'lineWidth', 2)
plot(ROI_boundaries(2, [1, 3]), ROI_boundaries(2, [2, 4]), 'r', 'lineWidth', 2)
plot(ROI_boundaries(3, [1, 3]), ROI_boundaries(3, [2, 4]), 'r', 'lineWidth', 2)
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
fig_wrapup(4, [])


%quantification
%1. including arena edges
center_score_entire = mean(mean(density_mat_e.*center_roi_quant));
periph_score_entire = mean(mean(density_mat_e.*periph_roi_quant));
entire_score = periph_score_entire./center_score_entire;

center_score_off = mean(mean(density_map_off.*center_roi_quant));
periph_score_off = mean(mean(density_map_off.*periph_roi_quant));
off_score = periph_score_off./center_score_off;


%2. excluding arena edges
center_score_entire = mean(mean(density_mat_e.*center_roi_quant_er));
periph_score_entire = mean(mean(density_mat_e.*periph_roi_quant_er));
entire_score_er = periph_score_entire./center_score_entire;

center_score_off = mean(mean(density_map_off.*center_roi_quant_er));
periph_score_off = mean(mean(density_map_off.*periph_roi_quant_er));
off_score_er = periph_score_off./center_score_off;


%3. statistical testing, testing null that distributions are equal on both
%sides of boundaries or ratio = 1
%entire arena
[h1, p1] = ttest(saved_scores(:, 1), 1);
%excluding edges
[h2, p2] = ttest(saved_scores(:, 2), 1);




function [] = play_trajs(traj_mat, mask_mat, traj_win, pulse_dur, frame_time)

n_trajs = size(traj_mat, 1)./2;
n_frames = size(traj_mat, 2);
another = 1;

while another == 1
    curr_traj_n = round(rand(1, 1).*(n_trajs - 1) + 1);
    curr_row_n = ((curr_traj_n - 1).*2) + 1;
    curr_traj = traj_mat(curr_row_n:(curr_row_n + 1), :);
    figure(3)
    plot_big_fig(3);
    imagesc(mask_mat)
    hold on
    for frame_n = 1:n_frames
        if frame_n < traj_win
            curr_color = [1, 1, 1];
        elseif frame_n >= traj_win && frame_n < (traj_win + pulse_dur)
            curr_color = [1, 0, 0];
        elseif frame_n >= (traj_win + pulse_dur)
            curr_color = [0, 0, 0];
        else
        end
        
        plot(curr_traj(1, frame_n), curr_traj(2, frame_n), '.', 'Color', curr_color, 'markerSize', 8);
        pause(frame_time);
        
    end    
    another = input('Play another trajectory? 1 - Yes, 0 - No');
    close figure 3;
end


end

function [occ_map] = generate_occ_map(traj_mat, n_interp)
occ_map = zeros(1280, 1024);
for traj_n = 1:(size(traj_mat, 1)./2)
    curr_traj = traj_mat(((traj_n - 1).*2 + 1):((traj_n - 1).*2 + 2), :);
    for frame_n = 1:(size(curr_traj, 2) - 1)
        curr_coords = curr_traj(:, frame_n:(frame_n + 1));
        curr_coords = round([linspace(curr_coords(1, 1), curr_coords(1, 2), n_interp); linspace(curr_coords(2, 1), curr_coords(2, 2), n_interp)]);     %interpolating between current pair of coords
        if sum(sum(isnan(curr_coords))) > 0
            continue
        else
        end
        
        occ_map(curr_coords(1, 1), curr_coords(2, 1)) = occ_map(curr_coords(1, 1), curr_coords(2, 1) ) + 1;     %loop below takes over from point2, so point1 done beforehand
        for interp_pt = 2:n_interp
            %ensuring there isn't double-counting due to interpolation
            if sum(curr_coords(:, interp_pt) - curr_coords(:, (interp_pt - 1))) ~= 0
                occ_map(curr_coords(1, interp_pt), curr_coords(2, interp_pt)) = occ_map(curr_coords(1, interp_pt), curr_coords(2, interp_pt) ) + 1;
            else
            end
            
        end
    end    
end

end
