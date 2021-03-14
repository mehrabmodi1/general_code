clear all
close all

base_order_paths = [{'C:\Data\Data\Adithya_arena_data\MB296B_TimeAirGap\CS+_First\'};...
             {'C:\Data\Data\Adithya_arena_data\MB296B_TimeAirGap\CS-_First\'}];

gap_paths = [{'Control\'}; {'15s\'}];
pulse_times_all = [{[31, 60;  61, 90]}; {[31, 60; 76, 105]}];
%pulse_times_all = [{[31, 60;  31, 60]}; {[31, 60; 31, 60]}];

%manually set parameters
equilib_time = 4;  %4;     %in s, the time for 1 vol-replacement in the arena because arena volume = pi*(5^2)*.3 = 23.6 cm^3 and flow rate = 400mL/min = 6.7 mL/s
t_window_orig = [0, 15]; %[0, 2]          %in s, manually chosen analysis time window after odor transition valve switch
%r_cutoff = [15, 35]; %[5, 45]   %in mm, the range of distances from center outside which flies are discarded as being too close to the arena center (0 mm) or edge (50 mm).
r_cutoff = [15, 40];   %in mm, the range of distances from center outside which flies are discarded as being too close to the arena center (0 mm) or edge (50 mm).

t_window_orig = t_window_orig + equilib_time;

score_vecs_all = [];
upwind_dist_tseries_all = [];
for base_o_path_n = 1:2
    base_order_path = base_order_paths{base_o_path_n};
    
    traj_samps_all = [];
    for gap_path_n = 1:2
        pulse_times = pulse_times_all{gap_path_n};
        t_window = t_window_orig + pulse_times(2, 1);        %interesting time point is onset of pulse2
        %t_window = t_window_orig + pulse_times(1, 2);         %looking at end of pulse1 - different for 15s gap datasets
        
        
        gap_path = gap_paths{gap_path_n};
        curr_path_base = [base_order_path, gap_path];
        dir_list = dir(curr_path_base);
        dir_list(1:2) = [];
            
        %loop to cycle through each experiment dataset
        score_vecs_gap = [];
        traj_samps_gap = [];
        upwind_dist_tseries_gap = [];
        for dir_n = 1:size(dir_list, 1)
            curr_dir = dir_list(dir_n).name;
            curr_path = [curr_path_base, curr_dir, '\'];
            curr_cami = findstr(curr_dir, '_Cam') + 4;
            curr_cam = curr_dir(curr_cami);

            %reading in metadata
            track_calib = load([curr_path, 'calibration.mat']);
            track_calib = track_calib.calib;
            frame_time = 1./track_calib.FPS;        %in s

            %reading in tracked data
            track_path = [curr_path, 'movie_Test_cam_', curr_cam, '\movie_Test_cam_', curr_cam, '-track.mat'];
            track_mat = load(track_path);
            track_mat = track_mat.trk;
            traj_mat = track_mat.data(:, :, 1:2);
%             if dir_n == 6
%                 keyboard
%             else
%             end
            
            traj_mat(:, :, 1) = traj_mat(:, :, 1) - max(track_calib.centroids);   %subtracting x-offset to set arena center to 0
            traj_mat(:, :, 2) = traj_mat(:, :, 2) - min(track_calib.centroids);   %subtracting y-offset to set arena center to 0             
            traj_mat = traj_mat./track_calib.PPM;       %converting position readings from pixels to mm
            
            
            %computing various behavioral scores
            %1. computing upwind dist travelled in t_window
            [traj_samps, upwind_dists, traj_mat, upwind_dist_tseries] = compute_center_dists(traj_mat, frame_time, t_window, r_cutoff);
            mean_upwind_vel = upwind_dists./(t_window(2) - t_window(1));     %mean upwind velocity in mm/s

            %2. computing total dist travelled in t_window
            [tot_dists] = compute_tot_dists(traj_mat, frame_time, t_window);
            
            score_name = 'upwind travel (mm)';
            score_vec = upwind_dists;

%             score_name = 'total travel (mm)';
%             score_vec = tot_dists;
            
%             score_name = 'frac. upwind travel';
%             score_vec = upwind_dists./tot_dists;
            
            score_vecs_gap = [score_vecs_gap; score_vec];
            traj_samps_gap = pad_n_concatenate(traj_samps_gap, traj_samps, 1, nan);
            upwind_dist_tseries_gap = pad_n_concatenate(upwind_dist_tseries_gap, upwind_dist_tseries, 1, nan);
        end
        score_vecs_all = pad_n_concatenate(score_vecs_all, score_vecs_gap, 2, nan);
        traj_samps_all = pad_n_concatenate(traj_samps_all, traj_samps_gap, 4, nan);
        upwind_dist_tseries_all = pad_n_concatenate(upwind_dist_tseries_all, upwind_dist_tseries_gap, 3, nan);
    end
    
    
    
    %PLOTTING
    %plotting sampled trajectories for 0 and 15s gaps
    end_colors = [0, 0, 0; 1, 0.4, 0.4;];
    %0s gap
%     plot_traj_samps(squeeze(traj_samps_all(:, :, :, 1)), end_colors, 2);
%     title('0s gap')
%     
%     plot_traj_samps(squeeze(traj_samps_all(:, :, :, 2)), end_colors, 3);
%     title('15s gap')
    
    
    keyboard
%     close figure 2
%     close figure 3
    
end

%plotting and statistical testing

%plotting distance time series' for flies
figure(4)
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 1));    %upwind distance time series for each valid fly for 0s gap
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', [0.1059, 0.6196, 0.4667]}, 1);
hold on
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 3));    %upwind distance time series for each valid fly for 0s gap
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', [0.8510, 0.3725, 0.0078]}, 1);
hold off

title('0s gap');
ylabel('upwind travel (mm)');
set_xlabels_time(4, frame_time, 5);
fig_wrapup(4, []);
ax_vals = axis;
ax_vals(3) = -4;
ax_vals(4) = 4;
ax_vals(2) = 59;
axis(ax_vals);




%plotting distance time series' for flies
figure(5)
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 2));    %upwind distance time series for each valid fly for 0s gap
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', [0.1059, 0.6196, 0.4667]}, 1);
hold on
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 4));    %upwind distance time series for each valid fly for 0s gap
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', [0.8510, 0.3725, 0.0078]}, 1);
hold off
title('15s gap');
ylabel('upwind travel (mm)');
set_xlabels_time(5, frame_time, 5);
fig_wrapup(5, []);
ax_vals = axis;
ax_vals(2) = 59;
ax_vals(3) = -4;
ax_vals(4) = 4;
axis(ax_vals);


unpaired_color = [0.1059, 0.6196, 0.4667];
paired_color = [0.8510, 0.3725, 0.0078];

score_vecs_all_final = [score_vecs_all(:, 1), score_vecs_all(:, 3), score_vecs_all(:, 2), score_vecs_all(:, 4)];        %re-arranging to bring paired, unpaired together instead of 0 and 15

markercolor = [unpaired_color; paired_color; unpaired_color.*0.7; paired_color.*0.7];
xlabels = [{'0 s, unpaired'}, {'0 s, paired'}, {'15 s, unpaired'}, {'15 s, paired'}];
fig_h = scattered_dot_plot_ttest(score_vecs_all_final, 1, 8, 8, 8, markercolor, 1, [], [], xlabels, 1, [0, 0, 0], 2, 0.05);
ylabel(score_name);
fig_wrapup(fig_h, []);




%-----------
%worker functions


function [] = plot_traj_samps(traj_samps, end_colors, fig_n)
    figure(fig_n);
    n_frames = size(traj_samps, 2);
    if sum(abs(end_colors(1, :) - end_colors(2, :))) ~= 0
        color_gradient = [linspace(end_colors(1, 1), end_colors(2, 1), n_frames)', linspace(end_colors(1, 2), end_colors(2, 2), n_frames)', linspace(end_colors(1, 3), end_colors(2, 3), n_frames)'];
        plot_gradient = 1;
    elseif sum(abs(end_colors(1, :) - end_colors(2, :))) == 0
        plot_gradient = 0;
    else
    end
    
    if plot_gradient == 1
        hold on
        for frame_n = 2:size(traj_samps, 2)
            curr_color = color_gradient(frame_n, :);
            plot(squeeze(traj_samps(:, (frame_n-1):frame_n, 1)'), squeeze(traj_samps(:, (frame_n-1):frame_n, 2)'), 'Color', curr_color);       
        end
        hold off
    elseif plot_gradient == 0
        curr_color = end_colors(1, :);
        plot(squeeze(traj_samps(:, :, 1)'), squeeze(traj_samps(:, :, 2)'), 'Color', curr_color);
    end
    xlabel('distance from center (mm)');
    ylabel('distance from center (mm)');
    fig_wrapup(fig_n, [])
    axis square
end

function [traj_samps, upwind_dists, traj_mat, upwind_dists_tseries] = compute_center_dists(traj_mat, frame_time, t_window, r_cutoff)
    loc_0 = squeeze(traj_mat(:, round(t_window(1)./frame_time), :));     %location at beginning of t_window, all flies
    edge_fliesi = find(sqrt(sum(loc_0.^2, 2)) < r_cutoff(1) | sqrt(sum(loc_0.^2, 2)) > r_cutoff(2));
    %edge_fliesi = find(sqrt(sum(loc_0.^2, 2)) < r_cutoff(2));
    traj_mat(edge_fliesi, :, :) = [];               %excluding flies too close to edge or center
    
    
    loc_0 = squeeze(traj_mat(:, round(t_window(1)./frame_time), :));     %location at beginning of t_window, non-edge flies
    loc_t = squeeze(traj_mat(:, round(t_window(2)./frame_time), :));     %location at end of t_window
    
    upwind_dists = sqrt(sum(loc_t.^2, 2)) - sqrt(sum(loc_0.^2, 2));      %distance from 0 ie. dist travelled upwind during t_window
    
    upwind_dists_tseries = [];
    
    
    for frame_n = (round(t_window(1)./frame_time) + 1):1:round(t_window(2)./frame_time)
        if frame_n == (round(t_window(1)./frame_time) + 1)
            zero_pts = squeeze(traj_mat(:, (frame_n - 1), :));
            zero_dists = sqrt(sum(zero_pts.^2, 2));
        else
        end
        curr_pts = squeeze(traj_mat(:, frame_n, :));
        curr_dists = sqrt(sum(curr_pts.^2, 2)) - zero_dists;           %distance from 0 in current frame
        upwind_dists_tseries = [upwind_dists_tseries, curr_dists];
    end
    pre_win = t_window(1) - (t_window(2) - t_window(1));        %sampling trajectories prior to beginning of t_window
    traj_samps = traj_mat(:, round(pre_win./frame_time):round(t_window(2)./frame_time), :);
    
end


function [tot_dists] = compute_tot_dists(traj_mat, frame_time, t_window)
    
    tot_dists = zeros(size(traj_mat, 1), 1);
    %loop to walk through each frame in t_window
    for frame_n = (round(t_window(1)./frame_time) + 1):1:round(t_window(2)./frame_time)
        curr_pos = squeeze(traj_mat(:, (frame_n), :));
        
        prev_pos = squeeze(traj_mat(:, (frame_n - 1), :));
        try
            curr_dists = sqrt( (curr_pos(:, 1) - prev_pos(:, 1)).^2 + (curr_pos(:, 2) - prev_pos(:, 2)).^2 );   %dist between positions in last and current frames for all flies
        catch
            keyboard
        end
        tot_dists = tot_dists + curr_dists;
    end
end
