clear all
close all

base_order_paths = [{'C:\Data\Data\Raw_data\Adithya_airgap_expts_v220210610\30s_Test\CS+_First\'};...
                    {'C:\Data\Data\Raw_data\Adithya_airgap_expts_v220210610\30s_Test\CS-_First\'}];
                
vid_path = 'C:\Data\Data\Analysed_data\Analysis_results\air_gap_traj_an\vert_aligned\';

gap_paths = [{'Control\'}; {'25s_gap\'}];
pulse_times_all = [{[31, 60;  61, 90]}; {[31, 60; 86, 115]}];  %real pulse-times
%pulse_times_all = [{[31, 60;  11, 30]}; {[31, 60; 11, 30]}];  %analyzing pre pulse1 baseline
%pulse_times_all = [{[31, 60;  61, 90]}; {[31, 60; 61, 90]}];  %analyzing pulse1 off period

analysis_offset = 0;
%analysis_offset = -10; %for delta norm.dist 0s gap
%analysis_offset = -35; %for delta norm.dist 25s gap

plot_video = 0;
align_vids_vert = 1;

vel_cutoffs = [2, .5];  %cutoffs in mm/s and s to determine if a fly has stopped or is moving

%splitting data by paired odor, or pooling all data if assigned empty
split_string = [];
%split_string = 'BA+';
%split_string = 'PA+';



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

cat_vec_baseline = mean(cat_vec(1:10000, 2));     %computing PID signal baseline (sampling first 1s)
cat_vec(:, 2) = cat_vec(:, 2) - cat_vec_baseline;

cat_vec_nogap_baseline = mean(cat_vec_nogap(1:10000, 2));     %computing PID signal baseline (sampling first 1s)
cat_vec_nogap(:, 2) = cat_vec_nogap(:, 2) - cat_vec_baseline;

frame_time = max(cat_vec(:, 1))./size(cat_vec, 1);
pad_t = (0:frame_time:25)';
pad = zeros(round(25./frame_time), 1) + nan;
cat_vec_nogap(:, 1) = cat_vec_nogap(:, 1) + 25;
cat_vec_nogap = [[pad_t, pad]; cat_vec_nogap];

figure(12)
t_offset = 91.5;
plot((cat_vec(:, 1) - t_offset), cat_vec(:, 2), 'Color', [0.65, 0.65, 0.65]);
hold on
plot((cat_vec_nogap(:, 1) - t_offset), cat_vec_nogap(:, 2), 'k');
ylabel('PID signal (V)')
xlabel('time (s)')
ax_vals = axis;
ax_vals(1) = -80;
ax_vals(2) = 48;
axis(ax_vals);
plot([0, 0], [ax_vals(3), ax_vals(4)], 'r')
fig_wrapup(12, [], [25, 30], .6);

%concluded that odor half-peak time from valve-switch is 6.5 s. Odor rise
%begins 5s after valve switch and reaches plateau in about 5 more s (total
%of 10s after valve switch). The 6.5s half-peak time is also valid for 0
%air-gap transitions.


%manually set parameters
equilib_time = 6.5;  %6.5; in s, Set as the time from valve-switch to reach odor half-peak. Time for 1 vol-replacement in the arena is ~4s because arena volume = pi*(5^2)*.3 = 23.6 cm^3 and flow rate = 400mL/min = 6.7 mL/s
t_window_orig = [0, 4]; %[0, 4]          %in s, manually chosen analysis time window after odor transition valve switch
%r_cutoff = [15, 35]; %[5, 45]   %in mm, the range of distances from center outside which flies are discarded as being too close to the arena center (0 mm) or edge (50 mm).

if pulse_times_all{2}(2, 1) == 61
    r_cutoff = [0, 50];    %don't want to exclude any flies for off-response analysis 
else
    r_cutoff = [10, 40];   %in mm, the range of distances from center outside which flies are discarded as being too close to the arena center (0 mm) or edge (50 mm).
end

%extending analysis window for analyses with offsets
if analysis_offset ~= 0 
    t_window_orig = [0, abs(analysis_offset)];
else
end

ang_cutoff = 0;      %in degrees, the angles at either end of the range that are excluded from analysis 
edge_cutoff = 48;     %in mm, the distance from center to identify flies whose upwind deviation scores should be highlighted

t_window_orig = t_window_orig + equilib_time + analysis_offset;

score_vecs_all = [];
downwind_deviations_all = [];
abs_dists_all = [];
edge_flies_all = [];
upwind_dist_tseries_all = [];
radial_pos_tseries_all = [];
downwind_deviations_tseries_all = [];
xy_vels_tseries_all = [];
xy_vels_bin_tseries_all = [];
ststp_tseries_all = [];
traj_mat_exps_all = [];
for base_o_path_n = 1:2
    base_order_path = base_order_paths{base_o_path_n};
    
    traj_samps_all = [];
    for gap_path_n = 1:2
        %skipping 0s gap for 25s gap offset analysis
        if analysis_offset == -35 && gap_path_n == 1
            t_window_orig_orig = t_window_orig;
            t_window_orig = [0, 4];
        elseif analysis_offset == -35 && gap_path_n == 2
            t_window_orig = t_window_orig_orig;
        end
        
        pulse_times = pulse_times_all{gap_path_n};
        t_window = t_window_orig + pulse_times(2, 1);        %interesting time point is onset of pulse2
        
        %t_window = t_window_orig + pulse_times(1, 2);         %looking at end of pulse1 - different for 25 s gap datasets
        
        
        gap_path = gap_paths{gap_path_n};
        curr_path_base = [base_order_path, gap_path];
        dir_list = dir(curr_path_base);
        dir_list(1:2) = [];
            
        %loop to cycle through each experiment dataset
        score_vecs_gap = [];
        downwind_deviations_vec = [];
        abs_dists_vec = [];
        edge_flies_vec = [];
        traj_samps_gap = [];
        upwind_dist_tseries_gap = [];
        radial_pos_tseries_gap = [];
        downwind_deviations_tserieses = [];
        xy_vels_tserieses = [];
        xy_vels_bin_tserieses = [];
        ststp_tserieses = [];
        traj_mat_exps = [];
        for dir_n = 1:size(dir_list, 1)
            curr_dir = dir_list(dir_n).name;
            curr_path = [curr_path_base, curr_dir, '\'];
            curr_cami = findstr(curr_dir, '_Cam') + 4;
            curr_cam = curr_dir(curr_cami);
            
            %splitting data by which odor was CS+
            if isempty(findstr(curr_path, split_string)) ~= 1      %only analyzing BA+
                continue
            else
            end
                       
            %reading in metadata
            track_calib = load([curr_path, 'calibration.mat']);
            track_calib = track_calib.calib;
            frame_time = 1./track_calib.FPS;        %in s

            %reading in tracked data
            track_path = [curr_path, 'movie_Test_cam_', curr_cam, '\movie_Test_cam_', curr_cam, '-track.mat'];
            track_mat = load(track_path);
            track_mat = track_mat.trk;
            traj_mat = track_mat.data(:, :, 1:3);
            
            traj_mat(:, :, 1) = traj_mat(:, :, 1) - max(track_calib.centroids);   %subtracting x-offset to set arena center to 0
            traj_mat(:, :, 2) = (traj_mat(:, :, 2) - min(track_calib.centroids)).* - 1;   %subtracting y-offset to set arena center to 0             
            traj_mat(:, :, 1:2) = traj_mat(:, :, 1:2)./track_calib.PPM;       %converting position readings from pixels to mm
            traj_mat_orig = traj_mat;
            
            %computing various behavioral scores
            %1. computing upwind dist travelled in t_window
            [traj_samps, upwind_dists, traj_mat, upwind_dist_tseries, traj_mat_exp, radial_pos_tseries] = compute_center_dists(traj_mat, frame_time, t_window, r_cutoff);
            mean_upwind_vel = upwind_dists./(t_window(2) - t_window(1));     %mean upwind velocity in mm/s
            
            %2. computing total dist travelled in t_window
            [tot_dists] = compute_tot_dists(traj_mat(:, :, 1:2), frame_time, t_window);
            score_name = 'upwind travel (mm)';
            score_vec = upwind_dists;
            
            %3. re-mapping cartesian orientations to radial orientations
            [downwind_deviations, downwind_deviations_tseries, edge_flies] = compute_radial_orientations(traj_mat_orig, frame_time, t_window, ang_cutoff, edge_cutoff);
            
            %4. computing mean distance from center over time window
            [mean_abs_dists] = compute_abs_dists(traj_mat, frame_time, t_window);
            
            %5. computing velocities in xy space (not upwind) and start-stop events based on velcity and run-duration thresholds.
            [xy_vels_tseries, xy_vels_bin, ststp_tseries] = get_xydists_ststp_probs(traj_mat, frame_time, t_window, vel_cutoffs);
            
            
            score_vecs_gap = [score_vecs_gap; score_vec];
            downwind_deviations_vec = [downwind_deviations_vec; downwind_deviations];
            abs_dists_vec = [abs_dists_vec; mean_abs_dists];
            edge_flies_vec = [edge_flies_vec; edge_flies];
            if size(traj_samps, 1) ~= 0
                traj_samps_gap = pad_n_concatenate(traj_samps_gap, traj_samps, 1, nan);
                upwind_dist_tseries_gap = pad_n_concatenate(upwind_dist_tseries_gap, upwind_dist_tseries, 1, nan);
                radial_pos_tseries_gap = pad_n_concatenate(radial_pos_tseries_gap, radial_pos_tseries, 1, nan);
                downwind_deviations_tserieses = pad_n_concatenate(downwind_deviations_tserieses, downwind_deviations_tseries, 1, nan);
                xy_vels_tserieses = pad_n_concatenate(xy_vels_tserieses, xy_vels_tseries, 1, nan);
                xy_vels_bin_tserieses = pad_n_concatenate(xy_vels_bin_tserieses, xy_vels_bin, 1, nan);
                ststp_tserieses = pad_n_concatenate(ststp_tserieses, ststp_tseries, 1, nan);
                traj_mat_exps = pad_n_concatenate(traj_mat_exps, traj_mat_exp, 1, nan);
            else
            end
            
        end
        score_vecs_all = pad_n_concatenate(score_vecs_all, score_vecs_gap, 2, nan);
        downwind_deviations_all = pad_n_concatenate(downwind_deviations_all, downwind_deviations_vec, 2, nan);
        abs_dists_all = pad_n_concatenate(abs_dists_all, abs_dists_vec, 2, nan);
        edge_flies_all = pad_n_concatenate(edge_flies_all, edge_flies_vec, 2, nan);
        traj_samps_all = pad_n_concatenate(traj_samps_all, traj_samps_gap, 4, nan);
        upwind_dist_tseries_all = pad_n_concatenate(upwind_dist_tseries_all, upwind_dist_tseries_gap, 3, nan);
        radial_pos_tseries_all = pad_n_concatenate(radial_pos_tseries_all, radial_pos_tseries_gap, 3, nan);
        downwind_deviations_tseries_all = pad_n_concatenate(downwind_deviations_tseries_all, downwind_deviations_tserieses, 3, nan);
        xy_vels_tseries_all = pad_n_concatenate(xy_vels_tseries_all, xy_vels_tserieses, 3, nan);
        xy_vels_bin_tseries_all = pad_n_concatenate(xy_vels_bin_tseries_all, xy_vels_bin_tserieses, 3, nan);
        ststp_tseries_all = pad_n_concatenate(ststp_tseries_all, ststp_tserieses, 3, nan);
        traj_mat_exps_all = pad_n_concatenate(traj_mat_exps_all, traj_mat_exps, 4, nan);
    end
    
end

%plotting and statistical testing
paired_color = [0,136,55]./256;
unpaired_color = [166,219,160]./256;

%plotting distance time series' for flies
figure(4)
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 3));    %upwind distance time series for each valid fly for 0s gap
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
hold on
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 1));    %upwind distance time series for each valid fly for 0s gap
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);

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
figure(5)
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 4));    %upwind distance time series for each valid fly for 0s gap
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
hold on
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 2));    %upwind distance time series for each valid fly for 0s gap
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);

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
%Plotting single traces
%plotting distance time series' for flies
figure(6)
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 3));    %upwind distance time series for each valid fly for 0s gap
plot(curr_traces', 'Color', paired_color);
hold on
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 1));    %upwind distance time series for each valid fly for 0s gap
plot(curr_traces', 'Color', unpaired_color);
hold off

title('0s gap');
ylabel('upwind travel (mm)');
set_xlabels_time(6, frame_time, 10);
fig_wrapup(6, [], [25, 30], .6);
ax_vals = axis;
ax_vals(3) = -25;
ax_vals(4) = 25;
axis(ax_vals);


%plotting distance time series' for flies
figure(7)
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 4));    %upwind distance time series for each valid fly for 0s gap
plot(curr_traces', 'Color', paired_color);
hold on
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 2));    %upwind distance time series for each valid fly for 0s gap
plot(curr_traces', 'Color', unpaired_color);
hold off
title('25 s gap');
ylabel('upwind travel (mm)');
set_xlabels_time(7, frame_time, 10);
fig_wrapup(7, [], [25, 30], .6);
ax_vals = axis;
ax_vals(3) = -25;
ax_vals(4) = 25;
axis(ax_vals);

%-----
%plotting normalized radial distance time series
%0s gap data
figure(22)
curr_traces = squeeze(radial_pos_tseries_all(:, :, 3));    %upwind distance time series for each valid fly for 0s gap
curr_traces = (curr_traces./50).^2;                         %normalizing to arena area at current distance

%Only done when computing delta for norm. distance with previous 4s
if analysis_offset == -10
    mid_f = floor(size(curr_traces, 2)./2);
    mid_c = mid_f + 1;
    curr_traces = curr_traces(:, mid_c:end) - curr_traces(:, 1:mid_f);
else
end
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
hold on
curr_traces = squeeze(radial_pos_tseries_all(:, :, 1));    %upwind distance time series for each valid fly for 0s gap
curr_traces = (curr_traces./50).^2;
%Only done when computing delta for norm. distance with previous 4s
if analysis_offset <= -10
    mid_f = floor(size(curr_traces, 2)./2);
    mid_c = mid_f + 1;
    curr_traces = curr_traces(:, mid_c:end) - curr_traces(:, 1:mid_f);
else
end
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);

hold off
title('0s gap');
ylabel('norm. distance');
set_xlabels_time(22, frame_time, 10);
fig_wrapup(22, [], [25, 30], .6);
ax_vals = axis;
ax_vals(3) = 0;
ax_vals(4) = .75;
%axis(2) = 
axis(ax_vals);

%25s gap data
figure(23)
curr_traces = squeeze(radial_pos_tseries_all(:, :, 4));    %upwind distance time series for each valid fly for 0s gap
curr_traces = (curr_traces./50).^2;                         %normalizing to arena area at current distance
%Only done when computing delta for norm. distance with previous 4s
if analysis_offset <= -10
    mid_f = floor(size(curr_traces, 2)./2);
    mid_c = mid_f + 1;
    curr_traces = curr_traces(:, mid_c:end) - curr_traces(:, 1:mid_f);
else
end
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
hold on
curr_traces = squeeze(radial_pos_tseries_all(:, :, 2));    %upwind distance time series for each valid fly for 0s gap
curr_traces = (curr_traces./50).^2;
%Only done when computing delta for norm. distance with previous 4s
if analysis_offset == -10
    mid_f = floor(size(curr_traces, 2)./2);
    mid_c = mid_f + 1;
    curr_traces = curr_traces(:, mid_c:end) - curr_traces(:, 1:mid_f);
else
end
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);

hold off
title('25 s gap');
ylabel('norm. distance');
set_xlabels_time(23, frame_time, 10);
fig_wrapup(23, [], [25, 30], .6);
ax_vals = axis;
ax_vals(3) = 0;
ax_vals(4) = .75;
axis(ax_vals);

%-----
%plotting upwind walking speeds
figure(25)
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 3));    %upwind distance time series for each valid fly for 0s gap
%computing speed
curr_traces = diff(curr_traces, 1, 2)./frame_time;       %speed in mm/s
curr_traces = movmean(curr_traces, 10, 2);
vel_mat(:, 1) = mean(curr_traces, 2, 'omitnan');

mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
hold on
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 1));    %upwind distance time series for each valid fly for 0s gap
%computing speed
curr_traces = diff(curr_traces, 1, 2)./frame_time;       %speed in mm/s
curr_traces = movmean(curr_traces, 10, 2);
vel_mat(:, 2) = mean(curr_traces, 2, 'omitnan');
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);

hold off
title('0s gap');
ylabel('upwind speed (mm/s)');
set_xlabels_time(25, frame_time, 10);
fig_wrapup(25, [], [25, 30], .6);
ax_vals = axis;
ax_vals(3) = -8;
ax_vals(4) = 8;
%ax_vals(2) = 59;
axis(ax_vals);


%plotting upwind walking speeds
figure(26)
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 4));    %upwind distance time series for each valid fly for 0s gap
%computing speed
curr_traces = diff(curr_traces, 1, 2)./frame_time;       %speed in mm/s
curr_traces = movmean(curr_traces, 10, 2);
vel_mat = pad_n_concatenate(vel_mat, mean(curr_traces, 2, 'omitnan'), 2, nan);
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
hold on
curr_traces = squeeze(upwind_dist_tseries_all(:, :, 2));    %upwind distance time series for each valid fly for 0s gap
%computing speed
curr_traces = diff(curr_traces, 1, 2)./frame_time;       %speed in mm/s
curr_traces = movmean(curr_traces, 10, 2);
vel_mat = pad_n_concatenate(vel_mat, mean(curr_traces, 2, 'omitnan'), 2, nan);
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);

hold off
title('25 s gap');
ylabel('upwind speed (mm/s)');
set_xlabels_time(26, frame_time, 10);
fig_wrapup(26, [], [25, 30], .6);
ax_vals = axis;
%ax_vals(2) = 59;
ax_vals(3) = -8;
ax_vals(4) = 8;
axis(ax_vals);


%plotting and statistical testing for upwind velocities
%plotting mean downwind deviation
score_vecs_all_final = [vel_mat(:, 2), vel_mat(:, 1), vel_mat(:, 4), vel_mat(:, 3)];      %re-arranging to bring paired, unpaired together instead of 0 and 15
markercolor = [unpaired_color; paired_color; unpaired_color; paired_color];
xlabels = [{'0 s, unprd'}, {'0 s, prd'}, {'25 s, unprd'}, {'25 s, prd'}];
fig_h = scattered_dot_plot_ttest(score_vecs_all_final, 27, 2.5, 4, 4, markercolor, 1, [], [], xlabels, 1, [0, 0, 0], 2, 0.05, 0);
ylabel('upwind speed (mm/s)');
fig_wrapup(fig_h, [], [25, 30], .6);



%plotting walking speed tseries
vel_mat = [];
figure(30)
curr_traces = squeeze(xy_vels_tseries_all(:, :, 3));    %upwind distance time series for each valid fly for 0s gap
curr_traces = movmean(curr_traces, 10, 2);
vel_mat(:, 1) = mean(curr_traces, 2, 'omitnan');
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
hold on
curr_traces = squeeze(xy_vels_tseries_all(:, :, 1));    %upwind distance time series for each valid fly for 0s gap
curr_traces = movmean(curr_traces, 10, 2);
vel_mat(:, 2) = mean(curr_traces, 2, 'omitnan');
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);

hold off
title('0s gap');
ylabel('speed (mm/s)');
set_xlabels_time(30, frame_time, 10);
fig_wrapup(30, [], [25, 30], .6);
ax_vals = axis;
ax_vals(3) = 0;
ax_vals(4) = 10;
%ax_vals(2) = 59;
axis(ax_vals);


%plotting walking speeds
figure(31)
curr_traces = squeeze(xy_vels_tseries_all(:, :, 4));    %upwind distance time series for each valid fly for 0s gap
curr_traces = movmean(curr_traces, 10, 2);
vel_mat = pad_n_concatenate(vel_mat, mean(curr_traces, 2, 'omitnan'), 2, nan);
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
hold on
curr_traces = squeeze(xy_vels_tseries_all(:, :, 2));    %upwind distance time series for each valid fly for 0s gap
curr_traces = movmean(curr_traces, 10, 2);
vel_mat = pad_n_concatenate(vel_mat, mean(curr_traces, 2, 'omitnan'), 2, nan);
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);

hold off
title('25 s gap');
ylabel('speed (mm/s)');
set_xlabels_time(31, frame_time, 10);
fig_wrapup(31, [], [25, 30], .6);
ax_vals = axis;
%ax_vals(2) = 59;
ax_vals(3) = 0;
ax_vals(4) = 10;
axis(ax_vals);


%plotting and statistical testing for undirected, xy velocities
%plotting mean running speed
score_vecs_all_final = [vel_mat(:, 2), vel_mat(:, 1), vel_mat(:, 4), vel_mat(:, 3)];      %re-arranging to bring paired, unpaired together instead of 0 and 15
markercolor = [unpaired_color; paired_color; unpaired_color; paired_color];
xlabels = [{'0 s, unprd'}, {'0 s, prd'}, {'25 s, unprd'}, {'25 s, prd'}];
fig_h = scattered_dot_plot_ttest(score_vecs_all_final, 32, 2.5, 4, 4, markercolor, 1, [], [], xlabels, 1, [0, 0, 0], 2, 0.05, 0);
ylabel('speed (mm/s)');
fig_wrapup(fig_h, [], [25, 30], .6);

%plotting running speed distribution
speed_vec = reshape(xy_vels_tseries_all, 1, []);
figure(34)
a = hist(speed_vec, [0:11]);
plot([0:10], a(1:11), 'lineWidth', 2);
title('xy speed histogram')
ylabel('counts')
xlabel('speed (mm/s)')
fig_wrapup(34, [], [25, 30], .6);



%computing and plotting start and stop probabilities
ststp_tseries_all_orig = ststp_tseries_all;
ststp_tseries_all(ststp_tseries_all < 0) = 0;   %measuring only start probabilities first
ststp_tseries_all = mean(ststp_tseries_all, 1, 'omitnan');

prob_mat = [];
figure(32)
curr_trace = squeeze(ststp_tseries_all(:, :, 3));    %upwind distance time series for each valid fly for 0s gap
curr_trace = movmean(curr_trace, 10, 2);
plot(curr_trace, 'lineWidth', 2, 'Color', paired_color);
hold on
curr_trace = squeeze(ststp_tseries_all(:, :, 1));    %upwind distance time series for each valid fly for 0s gap
curr_trace = movmean(curr_trace, 10, 2);
plot(curr_trace, 'lineWidth', 2, 'Color', unpaired_color);

hold off
title('0s gap');
ylabel('prob. transition');
set_xlabels_time(32, frame_time, 10);
fig_wrapup(32, [], [25, 30], .6);
ax_vals = axis;
ax_vals(3) = 0;
ax_vals(4) = 1;
%ax_vals(2) = 59;
axis(ax_vals);


%plotting transition probabilities
figure(33)
curr_trace = squeeze(ststp_tseries_all(:, :, 4));    %upwind distance time series for each valid fly for 0s gap
curr_trace = movmean(curr_trace, 10, 2);
plot(curr_trace, 'lineWidth', 2, 'Color', paired_color);
hold on
curr_trace = squeeze(ststp_tseries_all(:, :, 2));    %upwind distance time series for each valid fly for 0s gap
curr_trace = movmean(curr_trace, 10, 2);
prob_mat = pad_n_concatenate(prob_mat, mean(curr_traces, 2, 'omitnan'), 2, nan);
plot(curr_trace, 'lineWidth', 2, 'Color', unpaired_color);

hold off
title('25 s gap');
ylabel('prob. transition');
set_xlabels_time(33, frame_time, 10);
fig_wrapup(33, [], [25, 30], .6);
ax_vals = axis;
%ax_vals(2) = 59;
ax_vals(3) = 0;
ax_vals(4) = 1;
axis(ax_vals);



%plotting stop probabilities as dashed lines
ststp_tseries_all = ststp_tseries_all_orig;
ststp_tseries_all(ststp_tseries_all > 0) = 0;   %measuring only stop probabilities
ststp_tseries_all = mean(abs(ststp_tseries_all), 1, 'omitnan');

prob_mat = [];
figure(32)
hold on
curr_trace = squeeze(ststp_tseries_all(:, :, 3));    %upwind distance time series for each valid fly for 0s gap
curr_trace = movmean(curr_trace, 10, 2);
plot(curr_trace, ':', 'lineWidth', 2, 'Color', paired_color);
hold on
curr_trace = squeeze(ststp_tseries_all(:, :, 1));    %upwind distance time series for each valid fly for 0s gap
curr_trace = movmean(curr_trace, 10, 2);
plot(curr_trace, ':', 'lineWidth', 2, 'Color', unpaired_color);

hold off
title('0s gap');
ylabel('prob transition');
set_xlabels_time(32, frame_time, 10);
fig_wrapup(32, [], [25, 30], .6);
ax_vals = axis;
ax_vals(3) = 0;
ax_vals(4) = 0.05;
%ax_vals(2) = 59;
axis(ax_vals);



%plotting stop probabilities
figure(33)
hold on
curr_trace = squeeze(ststp_tseries_all(:, :, 4));    %upwind distance time series for each valid fly for 0s gap
curr_trace = movmean(curr_trace, 10, 2);
plot(curr_trace, ':', 'lineWidth', 2, 'Color', paired_color);
hold on
curr_trace = squeeze(ststp_tseries_all(:, :, 2));    %upwind distance time series for each valid fly for 0s gap
curr_trace = movmean(curr_trace, 10, 2);
prob_mat = pad_n_concatenate(prob_mat, mean(curr_traces, 2, 'omitnan'), 2, nan);
plot(curr_trace, ':', 'lineWidth', 2, 'Color', unpaired_color);

hold off
title('25 s gap');
ylabel('prob. transition');
set_xlabels_time(33, frame_time, 10);
fig_wrapup(33, [], [25, 30], .6);
ax_vals = axis;
%ax_vals(2) = 59;
ax_vals(3) = 0;
ax_vals(4) = 0.05;
axis(ax_vals);


%plotting probability of run v/s time
%plotting run probability
figure(35)
curr_trace = mean(squeeze(xy_vels_bin_tseries_all(:, :, 3)), 1, 'omitnan');    %upwind distance time series for each valid fly for 0s gap
curr_trace = movmean(curr_trace, 10, 2);
plot(curr_trace, 'lineWidth', 2, 'Color', paired_color);
hold on
curr_trace = mean(squeeze(xy_vels_bin_tseries_all(:, :, 1)), 1, 'omitnan');    %upwind distance time series for each valid fly for 0s gap
curr_trace = movmean(curr_trace, 10, 2);
plot(curr_trace, 'lineWidth', 2, 'Color', unpaired_color);

hold off
title('0s gap');
ylabel('prob. running');
set_xlabels_time(35, frame_time, 10);
fig_wrapup(35, [], [25, 30], .6);
ax_vals = axis;
ax_vals(3) = 0;
ax_vals(4) = 1;
%ax_vals(2) = 59;
axis(ax_vals);


%plotting transition probabilities
figure(36)
curr_trace = mean(squeeze(xy_vels_bin_tseries_all(:, :, 4)), 1, 'omitnan');    %upwind distance time series for each valid fly for 0s gap
curr_trace = movmean(curr_trace, 10, 2);
plot(curr_trace, 'lineWidth', 2, 'Color', paired_color);
hold on
curr_trace = mean(squeeze(xy_vels_bin_tseries_all(:, :, 2)), 1, 'omitnan');    %upwind distance time series for each valid fly for 0s gap
curr_trace = movmean(curr_trace, 10, 2);
prob_mat = pad_n_concatenate(prob_mat, mean(curr_traces, 2, 'omitnan'), 2, nan);
plot(curr_trace, 'lineWidth', 2, 'Color', unpaired_color);

hold off
title('25 s gap');
ylabel('prob. running');
set_xlabels_time(36, frame_time, 10);
fig_wrapup(36, [], [25, 30], .6);
ax_vals = axis;
%ax_vals(2) = 59;
ax_vals(3) = 0;
ax_vals(4) = 1;
axis(ax_vals);

%---------------



%PLOTTING downwind devIATION TIME SERIES
%plotting distance time series' for flies
%remapping upwind angles to have high abs degree values
downwind_deviations_tseries_all(downwind_deviations_tseries_all > 0) = 180 - downwind_deviations_tseries_all(downwind_deviations_tseries_all > 0);
downwind_deviations_tseries_all(downwind_deviations_tseries_all < 0) = -180 - downwind_deviations_tseries_all(downwind_deviations_tseries_all < 0);
% downwind_deviations_tseries_all(downwind_deviations_tseries_all > 0) = cos(deg2rad(abs(downwind_deviations_tseries_all(downwind_deviations_tseries_all > 0))));
% downwind_deviations_tseries_all(downwind_deviations_tseries_all < 0) = cos(deg2rad(abs(downwind_deviations_tseries_all(downwind_deviations_tseries_all < 0))));

figure(8)

curr_traces = abs(squeeze(downwind_deviations_tseries_all(:, :, 3)));    %downwind deviations time series for each valid fly for 0s gap
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
hold on
curr_traces = abs(squeeze(downwind_deviations_tseries_all(:, :, 1)));    %downwind deviations time series for each valid fly for 0s gap
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);

hold off
title('0s gap');
ylabel('upwind orientation (degrees)');
set_xlabels_time(8, frame_time, 10);
fig_wrapup(8, [], [25, 30], .6);
ax_vals = axis;
ax_vals(3) = -1;
ax_vals(4) = 1;
%ax_vals(2) = 59;
axis(ax_vals);


figure(9)
curr_traces = abs(squeeze(downwind_deviations_tseries_all(:, :, 4)));    %downwind deviationsance time series for each valid fly for 0s gap
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
hold on
curr_traces = abs(squeeze(downwind_deviations_tseries_all(:, :, 2)));    %downwind deviationsance time series for each valid fly for 0s gap
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);

hold off
title('25 s gap');
ylabel('upwind orientation (degrees)');
set_xlabels_time(9, frame_time, 10);
fig_wrapup(9, [], [25, 30], .6);
ax_vals = axis;
%ax_vals(2) = 59;
ax_vals(3) = 0;
ax_vals(4) = 180;
axis(ax_vals);


%------
%Plotting angular velocity
figure(28)

curr_traces = abs(squeeze(downwind_deviations_tseries_all(:, :, 3)));    %downwind deviations time series for each valid fly for 0s gap
curr_traces = diff(curr_traces, 1, 2)./frame_time;       %speed in degrees/s
curr_traces = movmean(curr_traces, 10, 2);
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
hold on
curr_traces = abs(squeeze(downwind_deviations_tseries_all(:, :, 1)));    %downwind deviations time series for each valid fly for 0s gap
curr_traces = diff(curr_traces, 1, 2)./frame_time;       %speed in degrees/s
curr_traces = movmean(curr_traces, 10, 2);
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);

hold off
title('0s gap');
ylabel('angular velocity (degrees/s)');
set_xlabels_time(28, frame_time, 10);
fig_wrapup(28, [], [25, 30], .6);
ax_vals = axis;
ax_vals(3) = -50;
ax_vals(4) = 50;
%ax_vals(2) = 59;
axis(ax_vals);


figure(29)
curr_traces = abs(squeeze(downwind_deviations_tseries_all(:, :, 4)));    %downwind deviationsance time series for each valid fly for 0s gap
curr_traces = diff(curr_traces, 1, 2)./frame_time;       %speed in degrees/s
curr_traces = movmean(curr_traces, 10, 2);
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
hold on
curr_traces = abs(squeeze(downwind_deviations_tseries_all(:, :, 2)));    %downwind deviationsance time series for each valid fly for 0s gap
curr_traces = diff(curr_traces, 1, 2)./frame_time;       %speed in degrees/s
curr_traces = movmean(curr_traces, 10, 2);
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);

hold off
title('25 s gap');
ylabel('angular velocity (degrees/s)');
set_xlabels_time(29, frame_time, 10);
fig_wrapup(29, [], [25, 30], .6);
ax_vals = axis;
%ax_vals(2) = 59;
ax_vals(3) = 0;
ax_vals(4) = 180;
axis(ax_vals);


%-------------



%plotting mean upwind displacement
score_vecs_all_final = [score_vecs_all(:, 1), score_vecs_all(:, 3), score_vecs_all(:, 2), score_vecs_all(:, 4)];        %re-arranging to bring paired, unpaired together instead of 0 and 15

markercolor = [unpaired_color; paired_color; unpaired_color; paired_color];
xlabels = [{'0 s, unprd'}, {'0 s, prd'}, {'25 s, unprd'}, {'25 s, prd'}];
fig_h = scattered_dot_plot_ttest(score_vecs_all_final, 1, 2.5, 4, 4, markercolor, 1, [], [], xlabels, 1, [0, 0, 0], 2, 0.05, 0);
ylabel(score_name);
fig_wrapup(fig_h, [], [25, 30], .6);

%plotting mean downwind deviation
score_vecs_all_final = 180 - [downwind_deviations_all(:, 1), downwind_deviations_all(:, 3), downwind_deviations_all(:, 2), downwind_deviations_all(:, 4)];        %re-arranging to bring paired, unpaired together instead of 0 and 15
markercolor = [unpaired_color; paired_color; unpaired_color; paired_color];
xlabels = [{'0 s, unprd'}, {'0 s, prd'}, {'25 s, unprd'}, {'25 s, prd'}];
fig_h = scattered_dot_plot_ttest(score_vecs_all_final, 2, 2.5, 4, 4, markercolor, 1, [], [], xlabels, 1, [0, 0, 0], 2, 0.05, 0);
ylabel('upwind orientation (degrees)');
fig_wrapup(fig_h, [], [25, 30], .6);


%plotting mean downwind deviation with edge flies highlighted
score_vecs_all_final = 180 - [downwind_deviations_all(:, 1), downwind_deviations_all(:, 3), downwind_deviations_all(:, 2), downwind_deviations_all(:, 4)];        %re-arranging to bring paired, unpaired together instead of 0 and 25
edge_vecs_all_final = [edge_flies_all(:, 1), edge_flies_all(:, 3), edge_flies_all(:, 2), edge_flies_all(:, 4)];        %re-arranging to bring paired, unpaired together instead of 0 and 25
markercolor = [unpaired_color; paired_color; unpaired_color; paired_color];
xlabels = [{'0 s, unprd'}, {'0 s, prd'}, {'25 s, unprd'}, {'25 s, prd'}];
[fig_h, r_vecs_saved] = scattered_dot_plot_ttest(score_vecs_all_final, 13, 2.5, 4, 4, markercolor, 1, [], [], xlabels, 1, [0, 0, 0], 2, 0.05, 0);
hold on
plot(r_vecs_saved(edge_vecs_all_final == 1), score_vecs_all_final(edge_vecs_all_final == 1), 'Or')
ylabel('upwind orientation (degrees)');
fig_wrapup(fig_h, [], [25, 30], .6);


%plotting mean absolute positions of flies
score_vecs_all_final = [abs_dists_all(:, 1), abs_dists_all(:, 3), abs_dists_all(:, 2), abs_dists_all(:, 4)];        %re-arranging to bring paired, unpaired together instead of 0 and 25
markercolor = [unpaired_color; paired_color; unpaired_color; paired_color];
xlabels = [{'0 s, unprd'}, {'0 s, prd'}, {'25 s, unprd'}, {'25 s, prd'}];
[fig_h, r_vecs_saved] = scattered_dot_plot_ttest(score_vecs_all_final, 14, 2.5, 4, 4, markercolor, 1, [], [], xlabels, 1, [0, 0, 0], 2, 0.05, 0);
hold on
ax_vals = axis;
plot([ax_vals(1), ax_vals(2)], [50, 50], 'r');
ylabel('dist. to center (mm)');
fig_wrapup(fig_h, [], [25, 30], .6);

%plotting mean norm. dists of flies
score_vecs_all_final = [abs_dists_all(:, 1), abs_dists_all(:, 3), abs_dists_all(:, 2), abs_dists_all(:, 4)];        %re-arranging to bring paired, unpaired together instead of 0 and 25
score_vecs_all_final = (score_vecs_all_final./50).^2;
markercolor = [unpaired_color; paired_color; unpaired_color; paired_color];
xlabels = [{'0 s, unprd'}, {'0 s, prd'}, {'25 s, unprd'}, {'25 s, prd'}];
[fig_h, r_vecs_saved] = scattered_dot_plot_ttest(score_vecs_all_final, 24, 2.5, 4, 4, markercolor, 1, [], [], xlabels, 1, [0, 0, 0], 2, 0.05, 0);
hold on
ax_vals = axis;
ylabel('norm. distance');
fig_wrapup(fig_h, [], [25, 30], .6);



%plotting starting position v/s upwind displacement for all flies and flies that pass r_cutoffs
base_path = 'C:\Data\Data\Analysed_data\Analysis_results\air_gap_traj_an\';
%reading in data
all_fly_displc = load([base_path, 'upwind_disps_allflies.mat']);
all_fly_displc = all_fly_displc.score_vecs_all_final;
sel_fly_displc = load([base_path, 'upwind_disps_selectflies.mat']);
sel_fly_displc = sel_fly_displc.score_vecs_all_final;

all_fly_pos = load([base_path, 'upwind_pos_allflies.mat']);
all_fly_pos = all_fly_pos.score_vecs_all_final;
sel_fly_pos = load([base_path, 'upwind_pos_selectflies.mat']);
sel_fly_pos = sel_fly_pos.score_vecs_all_final;

%0s gap
figure(37)
%unpaired odor as pulse2
plot(all_fly_pos(:, 1), all_fly_displc(:, 1), 'O', 'markerSize', 2.5, 'markerFaceColor', 'none', 'markerEdgeColor', unpaired_color);
hold on
plot(sel_fly_pos(:, 1), sel_fly_displc(:, 1), 'O', 'markerSize', 2.5, 'markerFaceColor', unpaired_color, 'markerEdgeColor', unpaired_color);

%paired odor as pulse2
plot(all_fly_pos(:, 3), all_fly_displc(:, 3), 'O', 'markerSize', 2.5, 'markerFaceColor', 'none', 'markerEdgeColor', paired_color);
plot(sel_fly_pos(:, 3), sel_fly_displc(:, 3), 'O', 'markerSize', 2.5, 'markerFaceColor', paired_color, 'markerEdgeColor', paired_color);
hold off

title('0s gap')
xlabel('start dist (mm)')
ylabel('upwind disp. (mm)')
ax_vals = axis;
ax_vals(1) = 0;
ax_vals(2) = 50;
ax_vals(3) = -40;
ax_vals(4) = 40;
plot([ax_vals(1), ax_vals(2)], [0, 0], 'Color', [0.65, 0.65, 0.65])
axis(ax_vals);

fig_wrapup(37, [], [25, 30], .6);


%25s gap
figure(38)
%unpaired odor as pulse2
plot(all_fly_pos(:, 2), all_fly_displc(:, 2), 'O', 'markerSize', 2.5, 'markerFaceColor', 'none', 'markerEdgeColor', unpaired_color);
hold on
plot(sel_fly_pos(:, 2), sel_fly_displc(:, 2), 'O', 'markerSize', 2.5, 'markerFaceColor', unpaired_color, 'markerEdgeColor', unpaired_color);

%paired odor as pulse2
plot(all_fly_pos(:, 4), all_fly_displc(:, 4), 'O', 'markerSize', 2.5, 'markerFaceColor', 'none', 'markerEdgeColor', paired_color);
plot(sel_fly_pos(:, 4), sel_fly_displc(:, 4), 'O', 'markerSize', 2.5, 'markerFaceColor', paired_color, 'markerEdgeColor', paired_color);
hold off
title('25s gap')
xlabel('start dist (mm)')
ylabel('upwind disp. (mm)')
ax_vals = axis;
ax_vals(1) = 0;
ax_vals(2) = 50;
ax_vals(3) = -40;
ax_vals(4) = 40;
axis(ax_vals);

fig_wrapup(38, [], [25, 30], .6);
    
    
 %off response analysis plots
if pulse_times_all{2}(2, 1) == 61      %analyzing off responses for 25s gap data
    
    %plotting upwind displacement for 25s gap data
    score_vecs_all_final = [score_vecs_all(:, 4), score_vecs_all(:, 2)];        %re-arranging to bring paired, unpaired together instead of 0 and 15
    markercolor = [unpaired_color; paired_color; unpaired_color; paired_color];
    xlabels = [{'25 s, unprd'}, {'25 s, prd'}];
    fig_h = scattered_dot_plot_ttest(score_vecs_all_final, 15, 2.5, 4, 4, markercolor, 1, [], [], xlabels, 1, [0, 0, 0], 2, 0.05, 0);
    title('pulse1 off response')
    ylabel(score_name);
    fig_wrapup(fig_h, [], [25, 30], .6);
    
    %plotting mean upwind orientation
    score_vecs_all_final = 180 - [downwind_deviations_all(:, 4), downwind_deviations_all(:, 3)];        %re-arranging to bring paired, unpaired together instead of 0 and 25
    edge_vecs_all_final = [edge_flies_all(:, 1), edge_flies_all(:, 3), edge_flies_all(:, 2), edge_flies_all(:, 4)];        %re-arranging to bring paired, unpaired together instead of 0 and 25
    markercolor = [unpaired_color; paired_color; unpaired_color; paired_color];
    xlabels = [{'25 s, unprd'}, {'25 s, prd'}];
    [fig_h, ~] = scattered_dot_plot_ttest(score_vecs_all_final, 16, 2.5, 4, 4, markercolor, 1, [], [], xlabels, 1, [0, 0, 0], 2, 0.05, 0);
    title('pulse1 off response')
    ylabel('upwind orientation (degrees)');
    fig_wrapup(fig_h, [], [25, 30], .6);

    
    %plotting absolute positions for 25s gap data   
    score_vecs_all_final = [abs_dists_all(:, 4), abs_dists_all(:, 2)];        %re-arranging to bring paired, unpaired together instead of 0 and 25
    markercolor = [unpaired_color; paired_color];
    xlabels = [{'25 s, unprd'}, {'25 s, prd'}];
    [fig_h, r_vecs_saved] = scattered_dot_plot_ttest(score_vecs_all_final, 17, 2.5, 4, 4, markercolor, 1, [], [], xlabels, 1, [0, 0, 0], 2, 0.05, 0);
    title('pulse1 off response')
    hold on
    ax_vals = axis;
    plot([ax_vals(1), ax_vals(2)], [50, 50], 'r');
    ylabel('dist. to center (mm)');
    fig_wrapup(fig_h, [], [25, 30], .6);
else
end


%plotting trajectory video
if plot_video == 1
    for dset_type = 1:4
        win_width = round( (t_window_orig(2) - t_window_orig(1))./frame_time);
        if dset_type == 1 | dset_type == 3  %0s gap case
            f0 = round( (t_window_orig(1) + pulse_times_all{1}(2, 1))./frame_time );
        elseif dset_type == 2 | dset_type == 4  %25s gap case
            f0 = round( (t_window_orig(1) + pulse_times_all{2}(2, 1))./frame_time );
        else
        end
        f1 = f0 - round(win_width./4);
        f2 = f0 + win_width;
        
        %computing time-stamp vector
        neg_vec = flipud((0:-frame_time:-((win_width./4).*frame_time))');    
        pos_vec = (0:frame_time:win_width.*frame_time)';                     
        t_vec = round([neg_vec; pos_vec(2:end)].*1000);     %in ms
        
        figure(17 + dset_type)
        set_names = [{'0s unprd'}, {'25s unprd'}, {'0s prd'}, {'25s prd'}];
        set_name = set_names(dset_type);
        set_name = set_name{1};
        v = VideoWriter([vid_path, set_name], 'MPEG-4');
        open(v)
        traj_mat = traj_mat_exps_all(:, f1:f2, :, dset_type);
        
        if align_vids_vert == 1
            %subtracting away t0 angular offset to align all flies to the upward vertical on frame1
            [thetaxy, rho] = cart2pol(traj_mat(:, :, 1), traj_mat(:, :, 2));
            thetaxy1 = thetaxy(:, 1);
            thetaxy = thetaxy - repmat(thetaxy1, 1, size(thetaxy, 2)) + pi./2;      %aligning frame1 angular positions to upper vertical
            [x, y] = pol2cart(thetaxy, rho);                        %going back to xy space
            
            %recomputing orientations to account for translation of thetaxy to upper vertical
            orientations = rad2deg(traj_mat(:, :, 3));
            del = find(orientations < 0);
            orientations(del) = 360 - abs(orientations(del));     %converting all degree measurements to positive values encoded as angles > 180 degrees
            orientations1 = orientations(:, 1);
            orientations = diff(orientations, 1, 2);              %expressing each orientation as a change from previous orientation
            orientations = [(orientations1 - rad2deg(thetaxy1)), orientations];    %adjusting orientation at t0 to account for reorientation at t0
            orientations = deg2rad(cumsum(orientations, 2)) + pi./2;                %adding orientation changes at each time point to newly established orientation at t0
            
            traj_mat(:, :, 3) = orientations;                               
            traj_mat(:, :, 1) = x;
            traj_mat(:, :, 2) = y;
        else
        end
        for frame_n = 1:(size(traj_mat, 2))
            for fly_n = 1:size(traj_mat, 1)
                curr_traj = traj_mat(fly_n, 1:end, :);
                
                
                
                if traj_mat_exps_all(fly_n, 1, 1, dset_type) == 0       %case where fly was excluded for being close to the edge/center
                    use_color = [0.9, 0.9, 0.9];
                elseif traj_mat_exps_all(fly_n, 1, 1, dset_type) == 1   %case where fly was included
                    use_color = [0, 0, 0];
                else
                end
                if isnan(curr_traj(1, frame_n, 1)) == 0
                    %computing curr_dist
                    curr_dist = mean(sqrt(sum(curr_traj(:, frame_n, 1:2).^2, 3)), 2, 'omitnan');
                    if curr_dist > 51
                        continue
                    else
                    end
                    
                    [pt2] = get_orientation_pt([curr_traj(1, frame_n, 1), curr_traj(1, frame_n, 2)], curr_traj(1, frame_n, 3), .5); %using tracked orientation to specify arrow direction
                    arr_pts = [curr_traj(1, (frame_n), 1), curr_traj(1, (frame_n), 2), pt2(1), pt2(2)];
                    plot_arrow(arr_pts(1), arr_pts(2), arr_pts(3), arr_pts(4), 'color', use_color, 'head_width', 2, 'head_height', 3, 'rel_sizing', 0, 'EdgeColor', [1, 1, 1], 'FaceColor', use_color);
                    t = text(35, -45, [num2str(t_vec(frame_n)), ' ms']);
                    t.FontSize = 12;
                    if t_vec(frame_n) >=0
                        t.Color = 'red';
                    else
                    end
                else
                end
                hold on

            end
            set(gca,'XColor', 'none','YColor','none');
            viscircles([0, 0], [51], 'lineWidth', 2, 'Color', [0, 0, 0]);
            viscircles([0, 0], [.1], 'lineWidth', 2, 'Color', [0, 0, 0]);
            axis equal
            axis square
            axis([-51, 51, -51, 51])
            set(gca,'XColor', 'none','YColor','none')
            
            drawnow
            hold off
            
            %writing fig frame to video file
            frame = getframe(gcf);
            writeVideo(v,frame);
           
        end
        keyboard
        close(v);
        
    end
    
    keyboard
    
else
end




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
    fig_wrapup(fig_n, [], [25, 30], .6)
    axis square
end

function [traj_samps, upwind_dists, traj_mat_ext, upwind_dists_tseries, traj_mat_exp, radial_pos_tseries] = compute_center_dists(traj_mat, frame_time, t_window, r_cutoff)
    
    %saving all trajectories, but logging which flies were not discarded as edge flies
    traj_mat_exp = zeros(size(traj_mat, 1), (size(traj_mat, 2) + 1), size(traj_mat, 3)) + 1;
    traj_mat_exp(:, 2:end, :) = traj_mat;
        
    traj_mat_ext = traj_mat;    
    traj_mat = traj_mat(:, :, 1:2);
    loc_0 = squeeze(traj_mat(:, round(t_window(1)./frame_time), :));     %location at beginning of t_window, all flies
    edge_fliesi = find(sqrt(sum(loc_0.^2, 2)) < r_cutoff(1) | sqrt(sum(loc_0.^2, 2)) > r_cutoff(2));
    %edge_fliesi = find(sqrt(sum(loc_0.^2, 2)) < r_cutoff(2));
    traj_mat(edge_fliesi, :, :) = [];               %excluding flies too close to edge or center
    traj_mat_ext(edge_fliesi, :, :) = [];
    traj_mat_exp(edge_fliesi, 1, :) = 0;
    
    loc_0 = squeeze(traj_mat(:, round(t_window(1)./frame_time), :));     %location at beginning of t_window, non-edge flies
    loc_t = squeeze(traj_mat(:, round(t_window(2)./frame_time), :));     %location at end of t_window
    
    upwind_dists = sqrt(sum(loc_t.^2, 2)) - sqrt(sum(loc_0.^2, 2));      %distance from 0 ie. dist travelled upwind during t_window
    
    upwind_dists_tseries = [];
    radial_pos_tseries = [];
    
    for frame_n = (round(t_window(1)./frame_time) + 1):1:round(t_window(2)./frame_time)
        if frame_n == (round(t_window(1)./frame_time) + 1)
            zero_pts = squeeze(traj_mat(:, (frame_n - 1), :));
            zero_dists = sqrt(sum(zero_pts.^2, 2));
        else
        end
        curr_pts = squeeze(traj_mat(:, frame_n, :));
        curr_dists_abs = sqrt(sum(curr_pts.^2, 2));
        curr_dists = curr_dists_abs - zero_dists;           %distance from 0 in current frame
        upwind_dists_tseries = [upwind_dists_tseries, curr_dists];
        radial_pos_tseries = [radial_pos_tseries, curr_dists_abs];
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

function [mean_downwind_deviations, downwind_deviations_tseries, edge_flies] = compute_radial_orientations(traj_mat, frame_time, t_window, ang_cutoff, edge_cutoff)
    
    radial_orientations = [];
    radial_orientations_tseries = [];

    %computing current position in arena in polar coordinates
    [theta, r] = cart2pol(traj_mat(:, :, 1), traj_mat(:, :, 2));
    theta = rad2deg(theta);
    ori_mat = rad2deg(squeeze(traj_mat(:, :, 3)));
    
    %correcting for angle convention flips around 180 degrees
    del = find(theta < 0);
    theta(del) = 360 - abs(theta(del));
    del = find(ori_mat < 0);
    ori_mat(del) = 360 - abs(ori_mat(del));
    
    %computing deviation from upwind orientation (0 = perfectly upwind orientation)
    radial_orientations_mat = ori_mat - theta;
    
    %correcting for abs(angles) > 180
    del = find(radial_orientations_mat > 180);
    radial_orientations_mat(del) = radial_orientations_mat(del) - 360;
    del = find(radial_orientations_mat < -180);
    radial_orientations_mat(del) = radial_orientations_mat(del) + 360;
    
    
%     %testing lines
%     fly_n = 15;
%     fr_n = [2000, 2700];
%         
%     figure(1)
%     plot(theta(fly_n, fr_n(1):fr_n(2))');
%     ylabel('pos angle');
%     
%     figure(2)
%     plot(ori_mat(fly_n, fr_n(1):fr_n(2))');
%     ylabel('tracked angle');
%     
%     figure(3)
%     plot(radial_orientations_mat(fly_n, fr_n(1):fr_n(2))');
%     ylabel('wind-rel angle');
%     
%     keyboard
       
    %sampling tseries in t_window
    fr1 = round(t_window(1)./frame_time);
    fr2 = round(t_window(2)./frame_time);
    
    downwind_deviations_tseries = radial_orientations_mat(:, fr1:fr2);
    
    %excluding flies whose orientations are too close to 0 or +/-180 at t = 0
%     ang_cutoffs = [0 - (ang_cutoff./2), 0 + (ang_cutoff./2)];       %for angles near 0
%     del = find(downwind_deviations_tseries(:, 1) > ang_cutoffs(1) & downwind_deviations_tseries(:, 1) < ang_cutoffs(2));
%     downwind_deviations_tseries(del, :) = []; 
%     traj_mat(del, :, :) = [];
%         
%     ang_cutoffs = [180 - (ang_cutoff./2), -180 + (ang_cutoff./2)];       %for angles near +/-180
%     del = find(downwind_deviations_tseries(:, 1) > ang_cutoffs(1) & downwind_deviations_tseries(:, 1) < ang_cutoffs(2));
%     downwind_deviations_tseries(del, :) = [];
%     traj_mat(del, :, :) = [];
    
    %identifying flies <2mm from upwind edge to highlight during plotting
    edge_flies = zeros(size(traj_mat, 1), 1);
    dists = mean(sqrt(sum(traj_mat(:, fr1:fr2, 1:2).^2, 3)), 2, 'omitnan');
    edge_flies(dists>edge_cutoff) = 1;
    
    %computing mean downind deviation
    mean_downwind_deviations = mean(abs(downwind_deviations_tseries), 2, 'omitnan');
    
end

function [mean_abs_dists] = compute_abs_dists(traj_mat, frame_time, t_window)
    %sampling tseries in t_window
    fr1 = round(t_window(1)./frame_time);
    fr2 = round(t_window(2)./frame_time);
    mean_abs_dists = mean(sqrt(sum(traj_mat(:, fr1:fr2, 1:2).^2, 3)), 2, 'omitnan');
end

function [pt2] = get_orientation_pt(pt1, angle, dist)
    %This function uses a point in xy space, and an orientation angle in
    %radians to specify a second point in xy space that is at dist and angle
    %relative to the first point.
    [x0, y0] = pol2cart(angle, dist);     %pt2 at specified angle and distance relative to the origin
    pt2 = [x0, y0] + pt1;

end

function [xy_vels_tseries, xy_vels_bin, ststp_tseries] = get_xydists_ststp_probs(traj_mat, frame_time, t_window, vel_cutoffs)

    t_window = round(t_window./frame_time);
    traj_mat = traj_mat(:, t_window(1):t_window(2), 1:2);
    
    %computing xy-velocities between consecutive pairs of points
    xy_vels_tseries = zeros(size(traj_mat, 1), size(traj_mat, 2));
    for fr_n = 2:size(traj_mat, 2)
        curr_pt1 = traj_mat(:, (fr_n - 1), :);
        curr_pt2 = traj_mat(:, fr_n, :);
        xy_vels_tseries(:, fr_n) = sum((curr_pt2 - curr_pt1).^2, 3).^0.5;     %xy distances between consecutive points in mm
    end
     xy_vels_tseries = xy_vels_tseries./frame_time;          %expressing instantaneous velocities in mm/s
    
    
     %Computing start/stop probabilities in time and run prob over time
     
    %binarizing xy-velocities around vel_cutoffs(1)
    xy_vels_bin = zeros(size(xy_vels_tseries, 1), size(xy_vels_tseries, 2));
    
    xy_vels_bin(xy_vels_tseries > vel_cutoffs(1)) = 1;                                  %running = 1, stopped = 0
    
    %getting rid of runs/stops that are shorter than vel_cutoffs(2)
    for fly_n = 1:size(xy_vels_bin, 1)
        xy_vels_bin(fly_n, :) = bwareaopen(xy_vels_bin(fly_n, :), round(vel_cutoffs(2)./frame_time));
    end    
    %inverting bin matrix
    xy_vels_bin = ones(size(xy_vels_bin, 1), size(xy_vels_bin, 2)) - xy_vels_bin;   %inverted, running = 0, stopped = 1
    %getting rid of stops that are too short
    for fly_n = 1:size(xy_vels_bin, 1)
        xy_vels_bin(fly_n, :) = bwareaopen(xy_vels_bin(fly_n, :), round(vel_cutoffs(2)./frame_time));
    end
    %re-inverting bin matrix to get running = 1, stopped = 0
    xy_vels_bin = ones(size(xy_vels_bin, 1), size(xy_vels_bin, 2)) - xy_vels_bin;   %binarised vel matrix, with runs/stops thresholded for a min duration
    
    
    %identifying start/stop transitions in binarised, thresholded vel matrix
    ststp_tseries = diff(xy_vels_bin, 1, 2);       %stp-run transitions = 1, run-stp ttransitions = -1
    
end

