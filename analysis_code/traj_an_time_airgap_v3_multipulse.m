clear all
close all

%multi pulse datasets
base_order_paths = [{'C:\Data\Data\Raw_data\Adithya_airgap_expts_v220210610\10sX6_Test\CS+_First\'};...
                    {'C:\Data\Data\Raw_data\Adithya_airgap_expts_v220210610\10sX6_Test\CS-_First\'}];

gap_paths = [{'Control\'}; {'25s_gap\'}];

%multi-pulse pulse-times
pulse_times_all = [{[31, 40; 41, 50; 51, 60; 61, 70; 71, 80; 81, 90; 91, 100; 101, 110; 111, 120; 121, 130; 131, 140; 141, 150]};...
                    {[31, 40; 66, 75; 101, 110; 136, 145; 171, 180; 206, 215; 241, 250; 276, 285; 311, 320; 346, 355; 381, 390; 416, 425]}];

use_pulse_set = [1, 2];
analysis_offset = 0;

%reading in PID traces to determine time adjustments to valve switch times.
pathgap = 'C:\Data\Data\Raw_data\Adithya_airgap_expts_v220210610\PID_Trace\CS+_AIR_CS-_30sTest\';
cat_vec = [];
for test_n = 1:3
    PID_trace = load([pathgap, 'PID_Test', num2str(test_n), '.mat']);
    PID_trace = PID_trace.PIDdata;
    if test_n > 1
        PID_trace(:, 1) = PID_trace(:, 1) + cat_vec(end, 1);
    else
    end
    cat_vec = [cat_vec; PID_trace(:, 1:2)];
end

pathnogap = 'C:\Data\Data\Raw_data\Adithya_airgap_expts_v220210610\PID_Trace\CS+_CS-_30sTest\';
cat_vec_nogap = [];
for test_n = 1:2
    PID_trace = load([pathnogap, 'PID_Test', num2str(test_n), '.mat']);
    PID_trace = PID_trace.PIDdata;
    if test_n > 1
        PID_trace(:, 1) = PID_trace(:, 1) + cat_vec_nogap(end, 1);
    else
    end
    cat_vec_nogap = [cat_vec_nogap; PID_trace(:, 1:2)];
end

figure(10)
plot(cat_vec(:, 1), cat_vec(:, 2));
hold on
plot(cat_vec_nogap(:, 1), cat_vec_nogap(:, 2), 'r');
ylabel('PID readout')
xlabel('time')

%concluded that odor half-peak time from valve-switch is 6.5 s. Odor rise
%begins 5s after valve switch and reaches plateau in about 5 more s (total
%of 10s after valve switch). The 6.5s half-peak time is also valid for 0
%air-gap transitions.

%manually set parameters
equilib_time = 6.5;  %6.5; in s, Set as the time from valve-switch to reach odor half-peak. Time for 1 vol-replacement in the arena is ~4s because arena volume = pi*(5^2)*.3 = 23.6 cm^3 and flow rate = 400mL/min = 6.7 mL/s
t_window_orig = [0, 4]; %[0, 2]          %in s, manually chosen analysis time window after odor transition valve switch
%r_cutoff = [15, 35]; %[5, 45]   %in mm, the range of distances from center outside which flies are discarded as being too close to the arena center (0 mm) or edge (50 mm).
r_cutoff = [10, 40];   %in mm, the range of distances from center outside which flies are discarded as being too close to the arena center (0 mm) or edge (50 mm).



t_window_orig = t_window_orig + equilib_time + analysis_offset;

score_vecs_all = [];
upwind_dist_tseries_all = [];
for base_o_path_n = 1:2
    base_order_path = base_order_paths{base_o_path_n};
    
    traj_samps_all = [];
    for gap_path_n = 1:2
        pulse_times_series = pulse_times_all{gap_path_n};
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
            
            score_vec1 = [];
                         
            for pulse_pair_n = use_pulse_set(1):2:use_pulse_set(2)
                pulse_times = pulse_times_series([pulse_pair_n, (pulse_pair_n + 1)], :);
                
                t_window = t_window_orig + pulse_times(2, 1);        %interesting time point is onset of pulse2
                %t_window = t_window_orig + pulse_times(1, 2);         %looking at end of pulse1 - different for 25 s gap datasets
                
                %computing various behavioral scores
                %1. computing upwind dist travelled in t_window
                
                [traj_samps, upwind_dists, traj_mat2, upwind_dist_tseries] = compute_center_dists(traj_mat, frame_time, t_window, r_cutoff);
                mean_upwind_vel = upwind_dists./(t_window(2) - t_window(1));     %mean upwind velocity in mm/s

                %2. computing total dist travelled in t_window
                [tot_dists] = compute_tot_dists(traj_mat2, frame_time, t_window);
                score_name = 'upwind travel (mm)';
                score_vec = upwind_dists;
                
                if isempty(score_vec) == 0
                    score_vecs_gap = [score_vecs_gap; score_vec];
                    traj_samps_gap = pad_n_concatenate(traj_samps_gap, traj_samps, 1, nan);
                    upwind_dist_tseries_gap = pad_n_concatenate(upwind_dist_tseries_gap, upwind_dist_tseries, 1, nan);
                else
                end
            end
        end
        score_vecs_all = pad_n_concatenate(score_vecs_all, score_vecs_gap, 2, nan);
        traj_samps_all = pad_n_concatenate(traj_samps_all, traj_samps_gap, 4, nan);
        upwind_dist_tseries_all = pad_n_concatenate(upwind_dist_tseries_all, upwind_dist_tseries_gap, 3, nan);
    end
    
    
    
    %PLOTTING
    %plotting sampled trajectories for 0 and 25 s gaps
    end_colors = [0, 0, 0; 1, 0.4, 0.4;];
    %0s gap
%     plot_traj_samps(squeeze(traj_samps_all(:, :, :, 1)), end_colors, 2);
%     title('0s gap')
%     
%     plot_traj_samps(squeeze(traj_samps_all(:, :, :, 2)), end_colors, 3);
%     title('25 s gap')
    
    
    %keyboard
%     close figure 2
%     close figure 3
    
end

%plotting and statistical testing
paired_color = [0,136,55]./256;
unpaired_color = [166,219,160]./256;

%plotting distance time series' for flies
figure(4)
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 1));    %upwind distance time series for each valid fly for 0s gap
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);
hold on
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 3));    %upwind distance time series for each valid fly for 0s gap
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
hold off

title('0s gap');
ylabel('upwind travel (mm)');
set_xlabels_time(4, frame_time, 5);
fig_wrapup(4, []);
ax_vals = axis;
ax_vals(3) = -4;
ax_vals(4) = 12;
%ax_vals(2) = 59;
axis(ax_vals);


%plotting distance time series' for flies
figure(5)
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 2));    %upwind distance time series for each valid fly for 0s gap
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);
hold on
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 4));    %upwind distance time series for each valid fly for 0s gap
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
hold off
title('25 s gap');
ylabel('upwind travel (mm)');
set_xlabels_time(5, frame_time, 5);
fig_wrapup(5, []);
ax_vals = axis;
%ax_vals(2) = 59;
ax_vals(3) = -4;
ax_vals(4) = 12;
axis(ax_vals);

%-----


%plotting distance time series' for flies
figure(6)
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 1));    %upwind distance time series for each valid fly for 0s gap
plot(curr_traces', 'Color', unpaired_color);
hold on
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 3));    %upwind distance time series for each valid fly for 0s gap
plot(curr_traces', 'Color', paired_color);
hold off

title('0s gap');
ylabel('upwind travel (mm)');
set_xlabels_time(4, frame_time, 10);
fig_wrapup(4, [], [25, 30], .6);
ax_vals = axis;
ax_vals(3) = -4;
ax_vals(4) = 12;
%ax_vals(2) = 59;
axis(ax_vals);


%plotting distance time series' for flies
figure(7)
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 2));    %upwind distance time series for each valid fly for 0s gap
plot(curr_traces', 'Color', unpaired_color);
hold on
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 4));    %upwind distance time series for each valid fly for 0s gap
plot(curr_traces', 'Color', paired_color);
hold off
title('25 s gap');
ylabel('upwind travel (mm)');
set_xlabels_time(5, frame_time, 10);
fig_wrapup(5, [], [25, 30], .6);
ax_vals = axis;
%ax_vals(2) = 59;
ax_vals(3) = -4;
ax_vals(4) = 12;
axis(ax_vals);

%-----


score_vecs_all_final = [score_vecs_all(:, 1), score_vecs_all(:, 3), score_vecs_all(:, 2), score_vecs_all(:, 4)];        %re-arranging to bring paired, unpaired together instead of 0 and 15

markercolor = [unpaired_color; paired_color; unpaired_color.*0.7; paired_color.*0.7];
xlabels = [{'0 s, unpaired'}, {'0 s, paired'}, {'25 s, unpaired'}, {'25 s, paired'}];
fig_h = scattered_dot_plot_ttest(score_vecs_all_final, 1, 2.5, 4, 4, markercolor, 1, [], [], xlabels, 1, [0, 0, 0], 2, 0.05, 0);
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
    if size(traj_mat, 1) == 1
        loc_0 = traj_mat(:, round(t_window(1)./frame_time), :);     %location at beginning of t_window, all flies
        loc_0 = reshape(loc_0, 1, size(loc_0, 3));     %getting rid of singleton dim2 (frame_n) but not singleton dim1 (fly_n)
    else
        loc_0 = squeeze(traj_mat(:, round(t_window(1)./frame_time), :));     %location at beginning of t_window, all flies
    end
        
    edge_fliesi = find(sqrt(sum(loc_0.^2, 2)) < r_cutoff(1) | sqrt(sum(loc_0.^2, 2)) > r_cutoff(2));
    %edge_fliesi = find(sqrt(sum(loc_0.^2, 2)) < r_cutoff(2));
    try
        traj_mat(edge_fliesi, :, :) = [];               %excluding flies too close to edge or center
    catch
        keyboard
    end
    
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
        
        %case when there is only 1 fly in traj_mat - squeeze screws up x-y
        %dim - need to flip
        if size(curr_pos, 2) == 1
            curr_pos = curr_pos';
            prev_pos = prev_pos';
        else
        end
        
        try
            curr_dists = sqrt( (curr_pos(:, 1) - prev_pos(:, 1)).^2 + (curr_pos(:, 2) - prev_pos(:, 2)).^2 );   %dist between positions in last and current frames for all flies
        catch
            keyboard
        end
        tot_dists = tot_dists + curr_dists;
    end
end
