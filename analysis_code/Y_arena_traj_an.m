clear all
close all

base_dir = 'C:\Data\Data\Adithya_Yarena_data\5sLED\';

dir_contents = dir(base_dir);
dir_contents(1:2) = [];
traj_win = 5;   %time window in s sampled for plotting
frame_time = 2./11;
traj_win = round(traj_win./frame_time);

traj_samples = [];
density_mat_e = zeros(1024, 1280);

%occupancy score matrix paramters
n_interp = 20;  %number of points to interpolate in between every pair of frames
blur_size = 2;  %SD of gaussian used to blur matrix

arena_ROI = load([base_dir, 'arena_ROI.mat']);
arena_ROI = arena_ROI.y;


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
    
    figure(1)
    set(gcf, 'Name', 'Entire trajectories overlaid');
    %getting stimulus defining ROIs
    masks = zeros(1024, 1280, 4);
    masks(:, :, 1) = data_struc.binarymask1;
    masks(:, :, 2) = data_struc.binarymask2 * 2;
    masks(:, :, 3) = data_struc.binarymask3 * 3;
    masks(:, :, 4) = data_struc.binarymask4 * 4;
    figure(1)
    if dir_n == 1
        imagesc(masks)
    else
    end
    
    hold on
    plot(coords(1, :), coords(2, :), '.w')
    hold on
%     led_on_pts = find(led_vec == 1);
%     plot(coords(1, led_on_pts), coords(2, led_on_pts), 'r.')
       
    %sampling traj segments around LED pulses
    d_led_vec = diff(led_vec);
    on_i = find(d_led_vec == 1);
    off_i = find(d_led_vec == -1);
    
    %accounting for interrupted LED pulse
    if length(off_i) == length(on_i) - 1
        on_i(end) = [];
    else
    end
    pulse_dur = median(off_i - on_i);
    
    for pulse_n = 1:size(on_i, 2)
        if (on_i(pulse_n) - traj_win) < 0 || (off_i(pulse_n) + traj_win) > size(coords, 2)
            continue
        else
        end
        traj_samp = coords(:, (on_i(pulse_n) - traj_win):(off_i(pulse_n) + traj_win));
        traj_samples = pad_n_concatenate(traj_samples, coords(:, (on_i(pulse_n) - traj_win):(off_i(pulse_n) + traj_win)), 1, nan);
        
    end
    
    %plotting trajectories sampled around LED pulse with post LED off
    %segments in black and pre LED on segments in white
    figure(2)
    imagesc(masks)
    hold on
%     plot(traj_samples((1:2:end), :), traj_samples((2:2:end), :), '.w')
%     plot(traj_samples((1:2:end), (traj_win + pulse_dur):end), traj_samples((2:2:end), (traj_win + pulse_dur):end), '.k')
%     plot(traj_samples((1:2:end), traj_win:(traj_win + pulse_dur) ), traj_samples((2:2:end), traj_win:(traj_win + pulse_dur) ), '.r')
%     
    
    %integrating into location density matrix for entire trajectory
    for frame_n = 1:(size(coords, 2) - 1)
        curr_coords = coords(:, frame_n:(frame_n + 1) );
        curr_coords = round([linspace(curr_coords(1, 1), curr_coords(1, 2), n_interp); linspace(curr_coords(2, 1), curr_coords(2, 2), n_interp)]);     %interpolating between current pair of coords
    
        for interp_pt = 1:n_interp
            density_mat_e(curr_coords(1, interp_pt), curr_coords(2, interp_pt)) = density_mat_e(curr_coords(1, interp_pt), curr_coords(2, interp_pt)) + 1;        
        end
    end
    

end

%play_trajs(traj_samples, masks, traj_win, pulse_dur, 0.075)

%integrating into density matrix for LED-sampled trajectories
density_mat_s = zeros(size(masks, 1), size(masks, 2));
for traj_n = 1:(size(traj_samples, 1)./2)
    curr_traj = traj_samples(((traj_n - 1).*2 + 1):((traj_n - 1).*2 + 2), :);
    for frame_n = 1:(size(curr_traj, 2) - 1)
        curr_coords = curr_traj(:, frame_n:(frame_n + 1));
        curr_coords = round([linspace(curr_coords(1, 1), curr_coords(1, 2), n_interp); linspace(curr_coords(2, 1), curr_coords(2, 2), n_interp)]);     %interpolating between current pair of coords
        if sum(sum(isnan(curr_coords))) > 0
            continue
        else
        end
        for interp_pt = 1:n_interp
            density_mat_s(curr_coords(1, interp_pt), curr_coords(2, interp_pt)) = density_mat_s(curr_coords(1, interp_pt), curr_coords(2, interp_pt) ) + 1;
            
        end
    end    
end


%making fly-averaged occupancy plot
figure(4)
set(gcf, 'Name', 'sampled trajectories overlaid')
imagesc(masks)
hold on
plot(traj_samples((1:2:end), :), traj_samples((2:2:end), :), '.w')
hold off


figure(5)
density_mat_e_bl = imgaussfilt(density_mat_e, blur_size);
set(gcf, 'Name', 'occupancy score plot, entire trajectories')
imagesc(density_mat_e_bl', [0, 4])

figure(6)
density_mat_s_bl = imgaussfilt(density_mat_s, blur_size);
set(gcf, 'Name', 'occupancy score plot, sampled trajectories')
imagesc(density_mat_s_bl', [0, 2])



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