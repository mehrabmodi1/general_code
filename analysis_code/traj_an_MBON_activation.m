clear all
close all


base_order_path = 'C:\Data\Data\Raw_data\Yoshi_MBON_activation\data\';
                    
                
vid_path = 'C:\Data\Data\Analysed_data\Analysis_results\air_gap_traj_an\vert_aligned\';
base_path = 'C:\Data\Data\Analysed_data\Analysis_results\air_gap_traj_an\';

pulse_width = 10;     %LED on duration in s

analysis_offset = -10;  %time before LED onset to start sampling trajectories
%analysis_offset = -10; %for delta norm.dist 0s gap
%analysis_offset = -35; %for delta norm.dist 25s gap

plot_video = 1;
align_vids_vert = 1;

included_trs = [5, 6];
average_trs = 1;
average_arenas = 1;     %0 - treat individual flies as samples; 1 - average flies within an arena and treat each arena as a sample
normalize_cent_dists = 0;      %(dist./arena radius).^2 ie. distances normalized by area at that distance
use_occupancy = 1;

vel_cutoffs = [2, .5];  %cutoffs in mm/s and s to determine if a fly has stopped or is moving
win_width = pulse_width + 18; %[0, 4]          %in s, manually chosen analysis time window width
%r_cutoff = [10, 40]; %[5, 45]   %in mm, the range of distances from center outside which flies are discarded as being too close to the arena center (0 mm) or edge (50 mm).
r_cutoff = [0, 50];        %computed inner cutoff such that the excluded inner area equals the excluded outer area of the arena


%manually set parameters
equilib_time = 0;  

ang_cutoff = 0;      %in degrees, the angles at either end of the range that are excluded from analysis 
edge_cutoff = 48;     %in mm, the distance from center to identify flies whose upwind deviation scores should be highlighted

score_vecs_all = [];
downwind_deviations_all = [];
abs_dists_all = [];
edge_flies_all = [];
upwind_dist_tseries_all = [];
radial_pos_tseries_all = [];
downwind_deviations_tseries_all = [];
downwind_deviations_tseries_singfly_all = [];
xydists_angles_tseries_all = [];
xy_vels_tseries_all = [];
xy_vels_bin_tseries_all = [];
ststp_tseries_all = [];
traj_mat_exps_all = [];


traj_samps_all = [];

curr_path_base = base_order_path;
dir_list = dir(curr_path_base);
dir_list(1:2) = [];

%loop to cycle through each experiment dataset
for dir_n = 1:size(dir_list, 1)
    curr_dir = dir_list(dir_n).name;
    curr_path = [curr_path_base, curr_dir, '\'];
    protocol = csvread([curr_path, 'protocol.csv'], 1, 0);
    %initial_delay = protocol(1, 20);
    initial_delay = 10;
    t_window = [initial_delay, (initial_delay + pulse_width)];
    
    upwind_dist_tseries_gap = [];
    radial_pos_tseries_gap = [];
    score_vecs_gap = [];
    downwind_deviations_vec = [];
    abs_dists_vec = [];
    edge_flies_vec = [];
    traj_samps_gap = [];


    downwind_deviations_tserieses = [];
    downwind_deviations_tserieses_singfly = [];
    xydists_angles_tserieses = [];
    xy_vels_tserieses = [];
    xy_vels_bin_tserieses = [];
    ststp_tserieses = [];
    traj_mat_exps = [];
    for LEDpulse_n = 1:6
        curr_cami = findstr(curr_dir, '_Cam') + 4;
        curr_cam = curr_dir(curr_cami);
        pulse_subpath = ['movie_AirLED', num2str(LEDpulse_n), '_cam_', num2str(curr_cam)];
        pulse_path = [curr_path, pulse_subpath, '\'];
        
        %reading in metadata
        track_calib = load([pulse_path, pulse_subpath, '-calibration.mat']);
        track_calib = track_calib.calib;
        frame_time = 1./track_calib.FPS;        %in s

        %reading in tracked data
        track_path = [pulse_path, pulse_subpath, '-track.mat'];
        track_mat = load(track_path);
        track_mat = track_mat.trk;
        traj_mat = track_mat.data(:, :, 1:3);

        traj_mat(:, :, 1) = traj_mat(:, :, 1) - max(track_calib.centroids);   %subtracting x-offset to set arena center to 0
        traj_mat(:, :, 2) = (traj_mat(:, :, 2) - min(track_calib.centroids)).* - 1;   %subtracting y-offset to set arena center to 0             
        traj_mat(:, :, 1:2) = traj_mat(:, :, 1:2)./track_calib.PPM;       %converting position readings from pixels to mm
        traj_mat_orig = traj_mat;

       
        %disp(['video duration = ' num2str(round(size(traj_mat, 2).*frame_time)), ' s']);
        %skipping the one extra-long video spotted by Yoshi
        if size(traj_mat, 2).*frame_time > 200
            continue
        else
        end

        %computing various behavioral scores
        %1. computing upwind dist travelled in t_window
        [traj_samps, upwind_dists, traj_mat, upwind_dist_tseries, traj_mat_exp, radial_pos_tseries] = compute_center_dists(traj_mat, frame_time, r_cutoff, normalize_cent_dists, equilib_time, pulse_width, win_width, initial_delay, use_occupancy);
        mean_upwind_vel = upwind_dists./pulse_width;     %mean upwind velocity in mm/s

        score_name = 'upwind travel (mm)';
        score_vec = upwind_dists;

        %3. re-mapping cartesian orientations to radial orientations
        t_window_2s = [initial_delay, (initial_delay + 2)];
        [downwind_deviations, downwind_deviations_tseries, edge_flies, xy_dists_angles_tseries] = compute_radial_orientations(traj_mat, frame_time, t_window_2s, ang_cutoff, edge_cutoff, analysis_offset);

        %4. computing mean distance from center over time window
        [mean_abs_dists] = compute_abs_dists(traj_mat, frame_time, t_window, analysis_offset);

        %5. computing velocities in xy space (not upwind) and start-stop events based on velcity and run-duration thresholds.
        [xy_vels_tseries, xy_vels_bin, ststp_tseries] = get_xydists_ststp_probs(traj_mat, frame_time, t_window, vel_cutoffs, analysis_offset);

        downwind_deviations_tseries_singfly = downwind_deviations_tseries;  %keeping a copy of measurements that need to be tracked for individual flies even when averaging within arenas

        %averaging within each arena if average_arenas == 1
        if average_arenas == 1
            score_vec = mean(score_vec, 1, 'omitnan');
            downwind_deviations = mean(downwind_deviations, 1, 'omitnan');
            mean_abs_dists = mean(mean_abs_dists, 1, 'omitnan');
            edge_flies = mean(edge_flies, 1, 'omitnan');

            if size(traj_samps, 1) ~= 0
                traj_samps = mean(traj_samps, 1, 'omitnan');
                upwind_dist_tseries = mean(upwind_dist_tseries, 1, 'omitnan');
                radial_pos_tseries = mean(radial_pos_tseries, 1, 'omitnan');
                downwind_deviations_tseries = mean(downwind_deviations_tseries, 1, 'omitnan');

                try
                    xy_vels_tseries = mean(xy_vels_tseries, 1, 'omitnan');

                catch
                    keyboard
                end
                xy_vels_bin = mean(xy_vels_bin, 1, 'omitnan');
                ststp_tseries = mean(ststp_tseries, 1, 'omitnan');
                traj_mat_exp = mean(traj_mat_exp, 1, 'omitnan');
            else
            end

        else
        end

        score_vecs_gap = [score_vecs_gap; score_vec];
        downwind_deviations_vec = [downwind_deviations_vec; downwind_deviations];
        abs_dists_vec = [abs_dists_vec; mean_abs_dists];
        edge_flies_vec = [edge_flies_vec; edge_flies];
        if size(traj_samps, 1) ~= 0
            traj_samps_gap = pad_n_concatenate(traj_samps_gap, traj_samps, 1, nan);
            upwind_dist_tseries_gap = pad_n_concatenate(upwind_dist_tseries_gap, upwind_dist_tseries, 1, nan);
            radial_pos_tseries_gap = pad_n_concatenate(radial_pos_tseries_gap, radial_pos_tseries, 1, nan);
            downwind_deviations_tserieses = pad_n_concatenate(downwind_deviations_tserieses, downwind_deviations_tseries, 1, nan);
            downwind_deviations_tserieses_singfly = pad_n_concatenate(downwind_deviations_tserieses_singfly, downwind_deviations_tseries_singfly, 1, nan);
            xydists_angles_tserieses = pad_n_concatenate(xydists_angles_tserieses, xy_dists_angles_tseries, 1, nan);
            try
                xy_vels_tserieses = pad_n_concatenate(xy_vels_tserieses, xy_vels_tseries, 1, nan);
            catch
                keyboard
            end
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
    downwind_deviations_tseries_singfly_all = pad_n_concatenate(downwind_deviations_tseries_singfly_all, downwind_deviations_tserieses_singfly, 3, nan);
    xydists_angles_tseries_all = pad_n_concatenate(xydists_angles_tseries_all, xydists_angles_tserieses, 3, nan);
    xy_vels_tseries_all = pad_n_concatenate(xy_vels_tseries_all, xy_vels_tserieses, 3, nan);
    xy_vels_bin_tseries_all = pad_n_concatenate(xy_vels_bin_tseries_all, xy_vels_bin_tserieses, 3, nan);
    ststp_tseries_all = pad_n_concatenate(ststp_tseries_all, ststp_tserieses, 3, nan);
    traj_mat_exps_all = pad_n_concatenate(traj_mat_exps_all, traj_mat_exps, 4, nan);
    
   
end


%Pooling the specified number of trials across datasets
upwind_dist_tseries_all_orig = upwind_dist_tseries_all;
subsamp = upwind_dist_tseries_all(included_trs(1):included_trs(2), :, :);
upwind_dist_tseries_all = [];
for dir_n = 1:size(subsamp, 3)
    if average_trs == 0
        curr_subsamp = squeeze(subsamp(:, :, dir_n));
    elseif average_trs == 1
        curr_subsamp = mean(squeeze(subsamp(:, :, dir_n)), 1, 'omitnan');
    end
    
    upwind_dist_tseries_all = [upwind_dist_tseries_all; curr_subsamp];    
end

%plotting and statistical testing
grey_col = [0.65, 0.65, 0.65];
red_col = [0.99, 0.06, 0.06];

%plotting distance time series' for flies
figure(2)
set(gcf, 'Name', 'mean upwind dists, 0s')
curr_traces = upwind_dist_tseries_all;     
mean_vec = mean(curr_traces, 1, 'omitnan');
se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
shadedErrorBar([], mean_vec, se_vec, {'Color', grey_col}, 1);
ylabel('upwind travel (mm)');
set_xlabels_time(2, frame_time, 10);
fig_wrapup(2, [], [25, 30], .6);
ax_vals = axis;
stim_frs = [round(10./frame_time), round(20./frame_time)];
fig_wrapup(2, [], [25, 30], 0.6);
add_stim_bar(2, stim_frs, red_col);

if normalize_cent_dists == 1
    y_level = 0.3;
else
    y_level = 7;
end
ax_vals(3) = -y_level;
ax_vals(4) = y_level;
%ax_vals(2) = 59;
axis(ax_vals);


%-----
%Plotting single traces
%plotting distance time series' for flies
figure(3)
set(gcf, 'Name', 'upwind dist traces, 0s')
curr_traces = squeeze(upwind_dist_tseries_all);     
plot(curr_traces', 'Color', grey_col);

ylabel('upwind travel (mm)');
set_xlabels_time(3, frame_time, 10);
fig_wrapup(3, [], [25, 30], .6);
ax_vals = axis;
if normalize_cent_dists == 1
    y_level = 0.5;
else
    y_level = 20;
end
add_stim_bar(3, stim_frs, red_col);
ax_vals(3) = -y_level;
ax_vals(4) = y_level;
axis(ax_vals);


%Pooling the specified number of trials across datasets
score_vecs_all_orig = score_vecs_all;
subsamp = score_vecs_all_orig(included_trs(1):included_trs(2), :);
score_vecs_all = [];
for dir_n = 1:size(subsamp, 2)
    if average_trs == 0
        curr_subsamp = subsamp(:, dir_n);
    elseif average_trs == 1
        curr_subsamp = mean(subsamp(:, dir_n), 1, 'omitnan');
    else
    end
    score_vecs_all = [score_vecs_all; curr_subsamp];    
end

%plotting mean upwind displacement
score_vecs_all_final = score_vecs_all;        
if average_trs == 0
    markercolor = [grey_col];
    xlabels = 'LED stim';
    fig_h = scattered_dot_plot_ttest(score_vecs_all_final, 4, .6, 1, 4, markercolor, 1, [], [], xlabels, 2, [0, 0, 0], 2, 0.05, 0, 1, 'force_mean');
elseif average_trs == 1
    figure(4)
    r_vec = zeros(size(score_vecs_all_final, 1)) + 1;
    fig_h = plot(r_vec, score_vecs_all_final, 'O', 'markerFaceColor', grey_col, 'markerEdgeColor', 'none', 'markerSize', 4);
    hold on
    mean_val = mean(score_vecs_all_final, 'omitnan');
    se_val = std(score_vecs_all_final, [], 'omitnan')./sqrt(size(score_vecs_all_final, 1));
    errorbar(1, mean_val, se_val, 'k', 'lineWidth', 2);
    plot(1, mean_val, 'O', 'markerSize', 5, 'markerFaceColor', [0, 0, 0], 'markerEdgeColor', 'none');
    hold off
else
end

ylabel(score_name);
set(gcf, 'Name', 'upwind displacement stats')
fig_wrapup(4, [], [25, 30], .6);
set(gca,'XTick',[])
[pval, hval] = signrank(score_vecs_all)

save('C:\Data\Data\Analysed_data\Analysis_results\air_gap_traj_an\MBON_LED_act_scores.mat', 'score_vecs_all');

%computing needed n to detect an effect size equal to LED activation mean, assuming normal distribution
%reference: https://www.itl.nist.gov/div898/handbook/prc/section2/prc222.htm
Z_val = 2.576;   %for alpha = 0.01, assuming normal distribution
SD = std(score_vecs_all, [], 'omitnan');
delta_detect = mean(score_vecs_all).*0.5;    %want to detect a difference in mean = (mean response with LED stim)./2
N = ((Z_val.*SD)./delta_detect).^2;

keyboard
% %fig_h = scattered_dot_plot_ttest(resp_mat_small, 4, 0.6, 1, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 1, 0.05, 0, 1, 'force_mean');
% 
% %Accounting for pulse1 off responses, but not odor-air responses
% %testing the null hypotheses that CS-on(0) == (CS+off(25) + CS-on(25)) and CS+on(0) == (CS-off(25) + CS+on(25))
% 
% %testing null hypotheses: Run these lines with pulse2 onset pulse times, as done usually
% pulse1_off_data = load([base_path, 'pulse1_off_upw_disps.mat']);
% pulse1_off_data = pulse1_off_data.score_vecs_all;
% 
% %bootstrapping 
% n_arenas = size(pulse1_off_data, 1);
% saved_null_pts = zeros(10000, 4);
% %1. Re-sampling n_arenas observations for each variable and computing null_pt_minus and null_pt_plus 10000 times 
% for r_samp_n = 1:10000
%     r_indices_off = round(rand(n_arenas, 4).*(n_arenas - 1)) + 1;   %sampling indices chosen with replacement
%     r_indices_on = round(rand(n_arenas, 4).*(n_arenas - 1)) + 1;   %sampling indices chosen with replacement
%     r_indices_curr = round(rand(n_arenas, 4).*(n_arenas - 1)) + 1;   %sampling indices chosen with replacement    
%     
%     r_pulse1_off_data = [pulse1_off_data(r_indices_off(:, 1), 1),  pulse1_off_data(r_indices_off(:, 2), 2)...
%                              pulse1_off_data(r_indices_off(:, 3), 3), pulse1_off_data(r_indices_off(:, 4), 4)];     %randomly re-sampling pulse1 off response data with replacement
%     r_pulse2_on_data = [score_vecs_all(r_indices_on(:, 1), 1),  score_vecs_all(r_indices_on(:, 2), 2)...
%                              score_vecs_all(r_indices_on(:, 3), 3), score_vecs_all(r_indices_on(:, 4), 4)];     %randomly re-sampling pulse1 off response data with replacement
%     r_pulse2_on_data_0s = [score_vecs_all(r_indices_curr(:, 1), 1),  score_vecs_all(r_indices_curr(:, 2), 2)...
%                              score_vecs_all(r_indices_curr(:, 3), 3), score_vecs_all(r_indices_curr(:, 4), 4)];     %randomly re-sampling pulse1 off response data with replacement
%         
%     CSplsoff25 = mean(r_pulse1_off_data(:, 4), 'omitnan');
%     CSmnsoff25 = mean(r_pulse1_off_data(:, 2), 'omitnan');
%     CSplson25 = mean(r_pulse2_on_data(:, 4), 'omitnan');
%     CSmnson25 = mean(r_pulse2_on_data(:, 2), 'omitnan');
% 
%     %plotting and statistical testing
%     null_pt_minus = CSplsoff25 + CSmnson25;     %null hypothesis test point CS-, pulse2 onset response with 0s gap
%     null_pt_plus = CSmnsoff25 + CSplson25;      %null hypothesis test point CS+, pulse2 onset response with 0s gap
%     
%     mean_responses_curr = mean(r_pulse2_on_data_0s(:, [1, 3]));     %mean 0s gap, upwind displacements
%     
%     mean_diffs = mean_responses_curr - [null_pt_minus, null_pt_plus];   %differences between actual means and linear model prediction means
%     
%     saved_null_pts(r_samp_n, :) = [null_pt_minus, null_pt_plus, mean_diffs];
%     
%     
% end
% null_pt_means = mean(saved_null_pts(:, 1:2), 'omitnan');
% null_pt_ses = std(saved_null_pts(:, 1:2), 'omitnan');       %STDs of re-sampled means approximate the SEMs of the null point distributions
% null_pt_var = (null_pt_ses.*sqrt(size(score_vecs_all, 1))).^2;     %STD = SE.*sqrt(n)
% curr_resp_means = mean(score_vecs_all(:, [1, 3]), 'omitnan');
% curr_resp_var = var(score_vecs_all(:, [1, 3]), 'omitnan')./sqrt(size(score_vecs_all, 1));
% 
% %computing two-sample t-statistic with means and ses of the two samples
% Sp_vals = sqrt((null_pt_var + curr_resp_var)./2);
% t_vals = (null_pt_means - curr_resp_means)./(Sp_vals.*sqrt(2./size(score_vecs_all, 1)));
% 
% p_vals = tcdf(t_vals, (size(score_vecs_all, 1) - 1));       %looking up p-values based on computed t-statistics
% 
% fig_h = scattered_dot_plot_ttest(score_vecs_all(:, [1,3]), 7, 2.5, 4, 4, markercolor, 1, [], [], xlabels, 2, [0, 0, 0], 2, 0.05, 0, 1, 'force_mean');
% 
% hold on
% 
% errorbar(1, null_pt_means(1), null_pt_ses(1), '.r', 'markerSize', 22, 'markerFaceColor', 'none');
% errorbar(2, null_pt_means(2), null_pt_ses(2), '.r', 'markerSize', 22, 'markerFaceColor', 'none');
% ylabel('upwind displacement (mm)');
% set(gcf, 'Name', 'upwind displacement stats')
% ax_vals = axis;
% ymax = 1.1.*ax_vals(4);
% text(0.9, ymax, ['p =', num2str(round(p_vals(1), 3))], 'FontSize', 7.5);
% text(1.9, ymax, ['p =',num2str(round(p_vals(2), 3))], 'FontSize', 7.5);
% fig_wrapup(fig_h, [], [25, 30], .6);
% 
% hold off
% keyboard
% %-----
% %plotting normalized radial distance time series
% %0s gap data
% % figure(7)
% % set(gcf, 'Name', 'Yoshi area norm, 0s')
% % curr_traces = squeeze(radial_pos_tseries_all(:, :, 3));     
% % curr_traces = (curr_traces./50).^2;                         %normalizing to arena area at current distance
% % 
% % %Only done when computing delta for norm. distance with previous 4s
% % if analysis_offset == -10
% %     mid_f = floor(size(curr_traces, 2)./2);
% %     mid_c = mid_f + 1;
% %     curr_traces = curr_traces(:, mid_c:end) - curr_traces(:, 1:mid_f);
% % else
% % end
% % mean_vec = mean(curr_traces, 1, 'omitnan');
% % se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
% % shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
% % hold on
% % curr_traces = squeeze(radial_pos_tseries_all(:, :, 1));     
% % curr_traces = (curr_traces./50).^2;
% % %Only done when computing delta for norm. distance with previous 4s
% % if analysis_offset <= -10
% %     mid_f = floor(size(curr_traces, 2)./2);
% %     mid_c = mid_f + 1;
% %     curr_traces = curr_traces(:, mid_c:end) - curr_traces(:, 1:mid_f);
% % else
% % end
% % mean_vec = mean(curr_traces, 1, 'omitnan');
% % se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
% % shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);
% % 
% % hold off
% % title('0s gap');
% % ylabel('norm. distance');
% % set_xlabels_time(6, frame_time, 10);
% % fig_wrapup(6, [], [25, 30], .6);
% % ax_vals = axis;
% % ax_vals(3) = 0;
% % ax_vals(4) = .75;
% % %axis(2) = 
% % axis(ax_vals);
% % 
% % %25s gap data
% % figure(7)
% % set(gcf, 'Name', 'Yoshi area norm, 25s')
% % 
% % curr_traces = squeeze(radial_pos_tseries_all(:, :, 4));     
% % curr_traces = (curr_traces./50).^2;                         %normalizing to arena area at current distance
% % %Only done when computing delta for norm. distance with previous 4s
% % if analysis_offset <= -10
% %     mid_f = floor(size(curr_traces, 2)./2);
% %     mid_c = mid_f + 1;
% %     curr_traces = curr_traces(:, mid_c:end) - curr_traces(:, 1:mid_f);
% % else
% % end
% % mean_vec = mean(curr_traces, 1, 'omitnan');
% % se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
% % shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
% % hold on
% % curr_traces = squeeze(radial_pos_tseries_all(:, :, 2));     
% % curr_traces = (curr_traces./50).^2;
% % %Only done when computing delta for norm. distance with previous 4s
% % if analysis_offset == -10
% %     mid_f = floor(size(curr_traces, 2)./2);
% %     mid_c = mid_f + 1;
% %     curr_traces = curr_traces(:, mid_c:end) - curr_traces(:, 1:mid_f);
% % else
% % end
% % mean_vec = mean(curr_traces, 1, 'omitnan');
% % se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
% % shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);
% % 
% % hold off
% % title('25 s gap');
% % ylabel('norm. distance');
% % set_xlabels_time(7, frame_time, 10);
% % fig_wrapup(7, [], [25, 30], .6);
% % ax_vals = axis;
% % ax_vals(3) = 0;
% % ax_vals(4) = .75;
% % axis(ax_vals);
% 
% %-----
% %plotting upwind walking speeds
% figure(8)
% set(gcf, 'Name', 'upwind vel, 0s')
% curr_traces = squeeze(upwind_dist_tseries_all(:, :, 3));     
% %computing speed
% curr_traces = diff(curr_traces, 1, 2)./frame_time;       %speed in mm/s
% curr_traces = movmean(curr_traces, 10, 2);
% vel_mat(:, 1) = mean(curr_traces, 2, 'omitnan');
% 
% mean_vec = mean(curr_traces, 1, 'omitnan');
% se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
% shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
% hold on
% curr_traces = squeeze(upwind_dist_tseries_all(:, :, 1));     
% %computing speed
% curr_traces = diff(curr_traces, 1, 2)./frame_time;       %speed in mm/s
% curr_traces = movmean(curr_traces, 10, 2);
% vel_mat(:, 2) = mean(curr_traces, 2, 'omitnan');
% mean_vec = mean(curr_traces, 1, 'omitnan');
% se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
% shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);
% 
% hold off
% title('0s gap');
% ylabel('upwind speed (mm/s)');
% set_xlabels_time(8, frame_time, 10);
% fig_wrapup(8, [], [25, 30], .6);
% ax_vals = axis;
% ax_vals(3) = -8;
% ax_vals(4) = 8;
% %ax_vals(2) = 59;
% axis(ax_vals);
% 
% 
% %plotting upwind walking speeds
% figure(9)
% set(gcf, 'Name', 'upwind vel, 25s')
% curr_traces = squeeze(upwind_dist_tseries_all(:, :, 4));     
% %computing speed
% curr_traces = diff(curr_traces, 1, 2)./frame_time;       %speed in mm/s
% curr_traces = movmean(curr_traces, 10, 2);
% vel_mat = pad_n_concatenate(vel_mat, mean(curr_traces, 2, 'omitnan'), 2, nan);
% mean_vec = mean(curr_traces, 1, 'omitnan');
% se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
% shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
% hold on
% curr_traces = squeeze(upwind_dist_tseries_all(:, :, 2));     
% %computing speed
% curr_traces = diff(curr_traces, 1, 2)./frame_time;       %speed in mm/s
% curr_traces = movmean(curr_traces, 10, 2);
% vel_mat = pad_n_concatenate(vel_mat, mean(curr_traces, 2, 'omitnan'), 2, nan);
% mean_vec = mean(curr_traces, 1, 'omitnan');
% se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
% shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);
% 
% hold off
% title('25 s gap');
% ylabel('upwind speed (mm/s)');
% set_xlabels_time(9, frame_time, 10);
% fig_wrapup(9, [], [25, 30], .6);
% ax_vals = axis;
% %ax_vals(2) = 59;
% ax_vals(3) = -8;
% ax_vals(4) = 8;
% axis(ax_vals);
% 
% 
% %plotting and statistical testing for upwind velocities
% %plotting mean downwind deviation
% score_vecs_all_final = [vel_mat(:, 2), vel_mat(:, 1), vel_mat(:, 4), vel_mat(:, 3)];      %re-arranging to bring paired, unpaired together instead of 0 and 15
% markercolor = [unpaired_color; paired_color; unpaired_color; paired_color];
% xlabels = [{'0 s, unprd'}, {'0 s, prd'}, {'25 s, unprd'}, {'25 s, prd'}];
% fig_h = scattered_dot_plot_ttest(score_vecs_all_final, 10, 2.5, 4, 4, markercolor, 1, [], [], xlabels, 2, [0, 0, 0], 2, 0.05, 0, 1, 'force_mean');
% set(gcf, 'Name', 'upwind vel stats, 0s')
% ylabel('upwind speed (mm/s)');
% fig_wrapup(fig_h, [], [25, 30], .6);
% 
% 
% xy_vels_tseries_all(xy_vels_tseries_all > 40) = nan;        %getting rid of junk values, ie, > 40 mm/s
% %plotting walking speed tseries
% vel_mat = [];
% figure(11)
% set(gcf, 'Name', 'xy walking speed, 0s')
% curr_traces = squeeze(xy_vels_tseries_all(:, :, 3));    %upwind vel time series for each valid fly for 0s gap
% curr_traces = movmean(curr_traces, 10, 2);
% vel_mat(:, 1) = mean(curr_traces, 2, 'omitnan');
% mean_vec = mean(curr_traces, 1, 'omitnan');
% se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
% shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
% hold on
% curr_traces = squeeze(xy_vels_tseries_all(:, :, 1));     
% curr_traces = movmean(curr_traces, 10, 2);
% vel_mat(:, 2) = mean(curr_traces, 2, 'omitnan');
% mean_vec = mean(curr_traces, 1, 'omitnan');
% se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
% shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);
% 
% hold off
% title('0s gap');
% ylabel('speed (mm/s)');
% set_xlabels_time(11, frame_time, 10);
% fig_wrapup(11, [], [25, 30], .6);
% ax_vals = axis;
% ax_vals(3) = 0;
% ax_vals(4) = 10;
% %ax_vals(2) = 59;
% axis(ax_vals);
% 
% 
% %plotting walking speeds
% figure(12)
% set(gcf, 'Name', 'xy walking speed, 25s')
% curr_traces = squeeze(xy_vels_tseries_all(:, :, 4));     
% curr_traces = movmean(curr_traces, 10, 2);
% vel_mat = pad_n_concatenate(vel_mat, mean(curr_traces, 2, 'omitnan'), 2, nan);
% mean_vec = mean(curr_traces, 1, 'omitnan');
% se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
% shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
% hold on
% curr_traces = squeeze(xy_vels_tseries_all(:, :, 2));     
% curr_traces = movmean(curr_traces, 10, 2);
% vel_mat = pad_n_concatenate(vel_mat, mean(curr_traces, 2, 'omitnan'), 2, nan);
% mean_vec = mean(curr_traces, 1, 'omitnan');
% se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
% shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);
% 
% hold off
% title('25 s gap');
% ylabel('speed (mm/s)');
% set_xlabels_time(11, frame_time, 10);
% fig_wrapup(12, [], [25, 30], .6);
% ax_vals = axis;
% %ax_vals(2) = 59;
% ax_vals(3) = 0;
% ax_vals(4) = 10;
% axis(ax_vals);
% 
% 
% %plotting and statistical testing for undirected, xy velocities
% %plotting mean running speed
% score_vecs_all_final = [vel_mat(:, 2), vel_mat(:, 1), vel_mat(:, 4), vel_mat(:, 3)];      %re-arranging to bring paired, unpaired together instead of 0 and 15
% markercolor = [unpaired_color; paired_color; unpaired_color; paired_color];
% xlabels = [{'0 s, unprd'}, {'0 s, prd'}, {'25 s, unprd'}, {'25 s, prd'}];
% fig_h = scattered_dot_plot_ttest(score_vecs_all_final, 13, 2.5, 4, 4, markercolor, 1, [], [], xlabels, 2, [0, 0, 0], 2, 0.05, 0, 1, 'force_mean');
% set(gcf, 'Name', 'xy vels stats, 0s')
% ylabel('speed (mm/s)');
% fig_wrapup(fig_h, [], [25, 30], .6);
% 
% %plotting running speed distribution
% % speed_vec = reshape(xy_vels_tseries_all, 1, []);
% % figure(14)
% % set(gcf, 'Name', 'xy speed hist')
% % a = hist(speed_vec, [0:11]);
% % plot([0:10], a(1:11), 'lineWidth', 2);
% % title('xy speed histogram')
% % ylabel('counts')
% % xlabel('speed (mm/s)')
% % fig_wrapup(14, [], [25, 30], .6);
% 
% 
% 
% %computing and plotting start and stop probabilities
% ststp_tseries_all_orig = ststp_tseries_all;
% ststp_tseries_all(ststp_tseries_all < 0) = 0;   %measuring only start probabilities first
% ststp_tseries_all = mean(ststp_tseries_all, 1, 'omitnan');
% 
% prob_mat = [];
% figure(15)
% set(gcf, 'Name', 'start probs, 0s')
% curr_trace = squeeze(ststp_tseries_all(:, :, 3));    
% curr_trace = movmean(curr_trace, 10, 2);
% plot(curr_trace, 'lineWidth', 2, 'Color', paired_color);
% hold on
% curr_trace = squeeze(ststp_tseries_all(:, :, 1));    
% curr_trace = movmean(curr_trace, 10, 2);
% plot(curr_trace, 'lineWidth', 2, 'Color', unpaired_color);
% 
% hold off
% title('0s gap');
% ylabel('prob. start');
% set_xlabels_time(15, frame_time, 10);
% fig_wrapup(15, [], [25, 30], .6);
% ax_vals = axis;
% ax_vals(3) = 0;
% ax_vals(4) = 0.05;
% %ax_vals(2) = 59;
% axis(ax_vals);
% 
% 
% %plotting transition probabilities
% figure(16)
% set(gcf, 'Name', 'start probs, 25s')
% curr_trace = squeeze(ststp_tseries_all(:, :, 4));     
% curr_trace = movmean(curr_trace, 10, 2);
% plot(curr_trace, 'lineWidth', 2, 'Color', paired_color);
% hold on
% curr_trace = squeeze(ststp_tseries_all(:, :, 2));     
% curr_trace = movmean(curr_trace, 10, 2);
% prob_mat = pad_n_concatenate(prob_mat, mean(curr_traces, 2, 'omitnan'), 2, nan);
% plot(curr_trace, 'lineWidth', 2, 'Color', unpaired_color);
% 
% hold off
% title('25 s gap');
% ylabel('prob. start');
% set_xlabels_time(16, frame_time, 10);
% fig_wrapup(16, [], [25, 30], .6);
% ax_vals = axis;
% %ax_vals(2) = 59;
% ax_vals(3) = 0;
% ax_vals(4) = 0.05;
% axis(ax_vals);
% 
% 
% 
% %plotting stop probabilities as dashed lines
% ststp_tseries_all = ststp_tseries_all_orig;
% ststp_tseries_all(ststp_tseries_all > 0) = 0;   %measuring only stop probabilities
% ststp_tseries_all = mean(abs(ststp_tseries_all), 1, 'omitnan');
% 
% prob_mat = [];
% figure(17)
% set(gcf, 'Name', 'stop probs, 0s')
% hold on
% curr_trace = squeeze(ststp_tseries_all(:, :, 3));     
% curr_trace = movmean(curr_trace, 10, 2);
% plot(curr_trace, 'lineWidth', 2, 'Color', paired_color);
% hold on
% curr_trace = squeeze(ststp_tseries_all(:, :, 1));     
% curr_trace = movmean(curr_trace, 10, 2);
% plot(curr_trace, 'lineWidth', 2, 'Color', unpaired_color);
% 
% hold off
% title('0s gap');
% ylabel('prob stop');
% set_xlabels_time(17, frame_time, 10);
% fig_wrapup(17, [], [25, 30], .6);
% ax_vals = axis;
% ax_vals(3) = 0;
% ax_vals(4) = 0.05;
% %ax_vals(2) = 59;
% axis(ax_vals);
% 
% 
% 
% %plotting stop probabilities
% figure(18)
% set(gcf, 'Name', 'stop probs, 25s')
% hold on
% curr_trace = squeeze(ststp_tseries_all(:, :, 4));     
% curr_trace = movmean(curr_trace, 10, 2);
% plot(curr_trace, 'lineWidth', 2, 'Color', paired_color);
% hold on
% curr_trace = squeeze(ststp_tseries_all(:, :, 2));     
% curr_trace = movmean(curr_trace, 10, 2);
% prob_mat = pad_n_concatenate(prob_mat, mean(curr_traces, 2, 'omitnan'), 2, nan);
% plot(curr_trace, 'lineWidth', 2, 'Color', unpaired_color);
% 
% hold off
% title('25 s gap');
% ylabel('prob. stop');
% set_xlabels_time(18, frame_time, 10);
% fig_wrapup(18, [], [25, 30], .6);
% ax_vals = axis;
% %ax_vals(2) = 59;
% ax_vals(3) = 0;
% ax_vals(4) = 0.05;
% axis(ax_vals);
% 
% 
% %plotting probability of run v/s time
% %plotting run probability
% figure(19)
% set(gcf, 'Name', 'prob. running, 0s')
% curr_trace = mean(squeeze(xy_vels_bin_tseries_all(:, :, 3)), 1, 'omitnan');    
% curr_trace = movmean(curr_trace, 10, 2);
% plot(curr_trace, 'lineWidth', 2, 'Color', paired_color);
% hold on
% curr_trace = mean(squeeze(xy_vels_bin_tseries_all(:, :, 1)), 1, 'omitnan');    
% curr_trace = movmean(curr_trace, 10, 2);
% plot(curr_trace, 'lineWidth', 2, 'Color', unpaired_color);
% 
% hold off
% title('0s gap');
% ylabel('prob. running');
% set_xlabels_time(19, frame_time, 10);
% fig_wrapup(19, [], [25, 30], .6);
% ax_vals = axis;
% ax_vals(3) = 0;
% ax_vals(4) = 1;
% %ax_vals(2) = 59;
% axis(ax_vals);
% 
% 
% %plotting transition probabilities
% figure(20)
% set(gcf, 'Name', 'prob. running, 25s')
% curr_trace = mean(squeeze(xy_vels_bin_tseries_all(:, :, 4)), 1, 'omitnan');     
% curr_trace = movmean(curr_trace, 10, 2);
% plot(curr_trace, 'lineWidth', 2, 'Color', paired_color);
% hold on
% curr_trace = mean(squeeze(xy_vels_bin_tseries_all(:, :, 2)), 1, 'omitnan');     
% curr_trace = movmean(curr_trace, 10, 2);
% prob_mat = pad_n_concatenate(prob_mat, mean(curr_traces, 2, 'omitnan'), 2, nan);
% plot(curr_trace, 'lineWidth', 2, 'Color', unpaired_color);
% 
% hold off
% title('25 s gap');
% ylabel('prob. running');
% set_xlabels_time(20, frame_time, 10);
% fig_wrapup(20, [], [25, 30], .6);
% ax_vals = axis;
% %ax_vals(2) = 59;
% ax_vals(3) = 0;
% ax_vals(4) = 1;
% axis(ax_vals);
% 
% %---------------
% 
% %PLOTTING downwind devIATION TIME SERIES
% %plotting distance time series' for flies
% %remapping upwind angles to have high abs degree values
% downwind_deviations_tseries_all(downwind_deviations_tseries_all > 0) = 180 - downwind_deviations_tseries_all(downwind_deviations_tseries_all > 0);
% downwind_deviations_tseries_all(downwind_deviations_tseries_all < 0) = -180 - downwind_deviations_tseries_all(downwind_deviations_tseries_all < 0);
% downwind_deviations_tseries_singfly_all(downwind_deviations_tseries_singfly_all > 0) = 180 - downwind_deviations_tseries_singfly_all(downwind_deviations_tseries_singfly_all > 0);
% downwind_deviations_tseries_singfly_all(downwind_deviations_tseries_singfly_all < 0) = -180 - downwind_deviations_tseries_singfly_all(downwind_deviations_tseries_singfly_all < 0);
% % downwind_deviations_tseries_all(downwind_deviations_tseries_all > 0) = cos(deg2rad(abs(downwind_deviations_tseries_all(downwind_deviations_tseries_all > 0))));
% % downwind_deviations_tseries_all(downwind_deviations_tseries_all < 0) = cos(deg2rad(abs(downwind_deviations_tseries_all(downwind_deviations_tseries_all < 0))));
% 
% figure(21)
% set(gcf, 'Name', 'upwind orientation, 0s')
% curr_traces = abs(squeeze(downwind_deviations_tseries_all(:, :, 3)));    %downwind deviations time series for each valid fly for 0s gap
% mean_vec = mean(curr_traces, 1, 'omitnan');
% se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
% shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
% hold on
% curr_traces = abs(squeeze(downwind_deviations_tseries_all(:, :, 1)));    %downwind deviations time series for each valid fly for 0s gap
% mean_vec = mean(curr_traces, 1, 'omitnan');
% se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
% shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);
% 
% hold off
% title('0s gap');
% ylabel('upwind orientation (degrees)');
% set_xlabels_time(21, frame_time, 10);
% fig_wrapup(21, [], [25, 30], .6);
% ax_vals = axis;
% ax_vals(3) = 0;
% ax_vals(4) = 180;
% %ax_vals(2) = 59;
% axis(ax_vals);
% 
% 
% figure(22)
% set(gcf, 'Name', 'upwind orientation, 25s')
% curr_traces = abs(squeeze(downwind_deviations_tseries_all(:, :, 4)));    %downwind deviationsance time series for each valid fly for 0s gap
% mean_vec = mean(curr_traces, 1, 'omitnan');
% se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
% shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
% hold on
% curr_traces = abs(squeeze(downwind_deviations_tseries_all(:, :, 2)));    %downwind deviationsance time series for each valid fly for 0s gap
% mean_vec = mean(curr_traces, 1, 'omitnan');
% se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
% shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);
% 
% hold off
% title('25 s gap');
% ylabel('upwind orientation (degrees)');
% set_xlabels_time(22, frame_time, 10);
% fig_wrapup(22, [], [25, 30], .6);
% ax_vals = axis;
% %ax_vals(2) = 59;
% ax_vals(3) = 0;
% ax_vals(4) = 180;
% axis(ax_vals);
% 
% 
% %plotting mean downwind deviation
% score_vecs_all_final = 180 - [downwind_deviations_all(:, 1), downwind_deviations_all(:, 3), downwind_deviations_all(:, 2), downwind_deviations_all(:, 4)];        %re-arranging to bring paired, unpaired together instead of 0 and 15
% markercolor = [unpaired_color; paired_color; unpaired_color; paired_color];
% xlabels = [{'0 s, unprd'}, {'0 s, prd'}, {'25 s, unprd'}, {'25 s, prd'}];
% fig_h = scattered_dot_plot_ttest(score_vecs_all_final, 23, 2.5, 4, 4, markercolor, 1, [], [], xlabels, 2, [0, 0, 0], 2, 0.05, 0, 1, 'force_mean');
% ylabel('upwind orientation (degrees)');
% set(gcf, 'Name', 'upwind orientation stats')
% fig_wrapup(fig_h, [], [25, 30], .6);
% 
% 
% %plotting mean downwind deviation with edge flies highlighted
% score_vecs_all_final = 180 - [downwind_deviations_all(:, 1), downwind_deviations_all(:, 3), downwind_deviations_all(:, 2), downwind_deviations_all(:, 4)];        %re-arranging to bring paired, unpaired together instead of 0 and 25
% edge_vecs_all_final = [edge_flies_all(:, 1), edge_flies_all(:, 3), edge_flies_all(:, 2), edge_flies_all(:, 4)];        %re-arranging to bring paired, unpaired together instead of 0 and 25
% markercolor = [unpaired_color; paired_color; unpaired_color; paired_color];
% xlabels = [{'0 s, unprd'}, {'0 s, prd'}, {'25 s, unprd'}, {'25 s, prd'}];
% [fig_h, r_vecs_saved] = scattered_dot_plot_ttest(score_vecs_all_final, 24, 2.5, 4, 4, markercolor, 1, [], [], xlabels, 2, [0, 0, 0], 2, 0.05, 0, 1, 'force_mean');
% hold on
% plot(r_vecs_saved(edge_vecs_all_final == 1), score_vecs_all_final(edge_vecs_all_final == 1), 'Or')
% ylabel('upwind orientation (degrees)');
% set(gcf, 'Name', 'upwind orientation, edge flies highlighted stats')
% fig_wrapup(fig_h, [], [25, 30], .6);
% 
% 
% %plotting mean downwind deviation
% score_vecs_all_final = 180 - [downwind_deviations_all(:, 1), downwind_deviations_all(:, 3), downwind_deviations_all(:, 2), downwind_deviations_all(:, 4)];        %re-arranging to bring paired, unpaired together instead of 0 and 15
% edgei = find(edge_vecs_all_final == 1);
% score_vecs_all_final(edgei) = nan;
% markercolor = [unpaired_color; paired_color; unpaired_color; paired_color];
% xlabels = [{'0 s, unprd'}, {'0 s, prd'}, {'25 s, unprd'}, {'25 s, prd'}];
% fig_h = scattered_dot_plot_ttest(score_vecs_all_final, 25, 2.5, 4, 4, markercolor, 1, [], [], xlabels, 2, [0, 0, 0], 2, 0.05, 0, 1, 'force_mean');
% ylabel('upwind orientation (degrees)');
% set(gcf, 'Name', 'upwind orientation stats, no edge flies')
% fig_wrapup(fig_h, [], [25, 30], .6);
% 
% %------
% %Plotting angular velocity
% figure(26)
% set(gcf, 'Name', 'angular velocity, 0s')
% curr_traces = abs(squeeze(downwind_deviations_tseries_all(:, :, 3)));    %downwind deviations time series for each valid fly for 0s gap
% curr_traces = diff(curr_traces, 1, 2)./frame_time;       %speed in degrees/s
% curr_traces = movmean(curr_traces, 10, 2);
% mean_vec = mean(curr_traces, 1, 'omitnan');
% se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
% shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
% hold on
% curr_traces = abs(squeeze(downwind_deviations_tseries_all(:, :, 1)));    %downwind deviations time series for each valid fly for 0s gap
% curr_traces = diff(curr_traces, 1, 2)./frame_time;       %speed in degrees/s
% curr_traces = movmean(curr_traces, 10, 2);
% mean_vec = mean(curr_traces, 1, 'omitnan');
% se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
% shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);
% 
% hold off
% title('0s gap');
% ylabel('angular velocity (degrees/s)');
% set_xlabels_time(26, frame_time, 10);
% fig_wrapup(26, [], [25, 30], .6);
% ax_vals = axis;
% ax_vals(3) = -50;
% ax_vals(4) = 50;
% %ax_vals(2) = 59;
% axis(ax_vals);
% 
% 
% figure(27)
% set(gcf, 'Name', 'angular velocity, 25s')
% curr_traces = abs(squeeze(downwind_deviations_tseries_all(:, :, 4)));    %downwind deviationsance time series for each valid fly for 0s gap
% curr_traces = diff(curr_traces, 1, 2)./frame_time;       %speed in degrees/s
% curr_traces = movmean(curr_traces, 10, 2);
% mean_vec = mean(curr_traces, 1, 'omitnan');
% se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
% shadedErrorBar([], mean_vec, se_vec, {'Color', paired_color}, 1);
% hold on
% curr_traces = abs(squeeze(downwind_deviations_tseries_all(:, :, 2)));    %downwind deviationsance time series for each valid fly for 0s gap
% curr_traces = diff(curr_traces, 1, 2)./frame_time;       %speed in degrees/s
% curr_traces = movmean(curr_traces, 10, 2);
% mean_vec = mean(curr_traces, 1, 'omitnan');
% se_vec = std(curr_traces, [], 1, 'omitnan')./sqrt(size(curr_traces, 1));
% shadedErrorBar([], mean_vec, se_vec, {'Color', unpaired_color}, 1);
% 
% hold off
% title('25 s gap');
% ylabel('angular velocity (degrees/s)');
% set_xlabels_time(27, frame_time, 10);
% fig_wrapup(27, [], [25, 30], .6);
% ax_vals = axis;
% %ax_vals(2) = 59;
% ax_vals(3) = -50;
% ax_vals(4) = 50;
% axis(ax_vals);
% 
% 
% %-------------






%plotting mean absolute positions of flies
score_vecs_all_final = [abs_dists_all(:, 1), abs_dists_all(:, 3), abs_dists_all(:, 2), abs_dists_all(:, 4)];        %re-arranging to bring paired, unpaired together instead of 0 and 25
markercolor = [unpaired_color; paired_color; unpaired_color; paired_color];
xlabels = [{'0 s, unprd'}, {'0 s, prd'}, {'25 s, unprd'}, {'25 s, prd'}];
[fig_h, r_vecs_saved] = scattered_dot_plot_ttest(score_vecs_all_final, 28, 2.5, 4, 4, markercolor, 1, [], [], xlabels, 2, [0, 0, 0], 2, 0.05, 0, 1, 'force_mean');
hold on
ax_vals = axis;
plot([ax_vals(1), ax_vals(2)], [50, 50], 'r');
set(gcf, 'Name', 'mean dist to center stats')
ylabel('dist. to center (mm)');
fig_wrapup(fig_h, [], [25, 30], .6);

if (t_window_orig(2) - t_window_orig(1)) == 4
    save([base_path, 'upwind_disps_allflies.mat'], 'abs_dists_all');
elseif (t_window_orig(2) - t_window_orig(1)) < 0.21
    save([base_path, 'upwind_pos_allflies.mat'], 'abs_dists_all');
else
end

% %Plotting orientation v/s xy speed for all time-points in analysis window
% figure(35)
% set(gcf, 'Name', 'downwind deviations v/s xy speeds, 0s')
% %grouping xyspeeds according to their corressponding orientation angles, for each fly
% angle_bins = 0:20:180;
% downwind_deviations_tseries_singfly_all = abs(downwind_deviations_tseries_singfly_all);
% all_speeds_ang = [];
% 
% %paired odor, 0s
% mean_spd_vecs_all = zeros(size(downwind_deviations_tseries_singfly_all, 1), (length(angle_bins)) - 1) + nan;
% %mean_angvel_vecs_all = zeros(size(downwind_deviations_tseries_all, 1), (length(angle_bins)) - 1) + nan;
% for fly_n = 1:size(downwind_deviations_tseries_singfly_all, 1)
%     curr_devs = downwind_deviations_tseries_singfly_all(fly_n, :, 1);
%     xy_speeds = xydists_angles_tseries_all(fly_n, :, 1)./frame_time;
%     [grouped_vals_vec2] = bin_vec2_by_vec1(curr_devs, xy_speeds, angle_bins);   %xyspeeds grouped by corress heading angle for a given fly
%     mean_spd_vecs_all(fly_n, :) = mean(grouped_vals_vec2, 1, 'omitnan');
% end
% 
% mean_spds = mean(mean_spd_vecs_all, 1, 'omitnan');                                          %average across flies
% se_spds = std(mean_spd_vecs_all, [], 1, 'omitnan')./sqrt(size(mean_spd_vecs_all, 1));       %SE across flies
% shadedErrorBar((angle_bins(2:end) - 10), mean_spds, se_spds, {'Color', paired_color}, 1);   %subtracting 10 from bins to set to middle rather than upper edge of each bin
% hold on
% 
% %logging grouped speeds for statistical testing
% speeds_ang = pad_n_concatenate(reshape(mean_spd_vecs_all(:, 1:3), [], 1), reshape(mean_spd_vecs_all(:, 7:9), [], 1), 2, nan);
% all_speeds_ang = pad_n_concatenate(all_speeds_ang, speeds_ang, 2, nan);
% 
% 
% %unpaired odor, 0s
% mean_spd_vecs_all = zeros(size(downwind_deviations_tseries_singfly_all, 1), (length(angle_bins)) - 1) + nan;
% for fly_n = 1:size(downwind_deviations_tseries_singfly_all, 1)
%     curr_devs = downwind_deviations_tseries_singfly_all(fly_n, :, 3);
%     xy_speeds = xydists_angles_tseries_all(fly_n, :, 3)./frame_time;
%     [grouped_vals_vec2] = bin_vec2_by_vec1(curr_devs, xy_speeds, angle_bins);   %xyspeeds grouped by corress heading angle for a given fly
%     mean_spd_vecs_all(fly_n, :) = mean(grouped_vals_vec2, 1, 'omitnan');
% end
% mean_spds = mean(mean_spd_vecs_all, 1, 'omitnan');
% se_spds = std(mean_spd_vecs_all, [], 1, 'omitnan')./sqrt(size(mean_spd_vecs_all, 1));
% shadedErrorBar((angle_bins(2:end) - 10), mean_spds, se_spds, {'Color', unpaired_color}, 1);
% hold off
% 
% %logging grouped speeds for statistical testing
% speeds_ang = pad_n_concatenate(reshape(mean_spd_vecs_all(:, 1:3), [], 1), reshape(mean_spd_vecs_all(:, 7:9), [], 1), 2, nan);
% all_speeds_ang = pad_n_concatenate(all_speeds_ang, speeds_ang, 2, nan);
% 
% 
% title('0s gap')
% xlabel('abs. deviation from upwind (degrees)')
% ylabel('xy speed (mm/s)')
% ax_vals = axis;
% ax_vals(1) = 0;
% ax_vals(2) = 180;
% ax_vals(3) = 0;
% ax_vals(4) = 15;
% axis(ax_vals);
% 
% fig_wrapup(35, [], [25, 30], .6);
% 
% figure(36)
% set(gcf, 'Name', 'downwind deviations v/s xy speeds, 25s')
% %paired odor, 25s
% mean_spd_vecs_all = zeros(size(downwind_deviations_tseries_singfly_all, 1), (length(angle_bins)) - 1) + nan;
% for fly_n = 1:size(downwind_deviations_tseries_singfly_all, 1)
%     curr_devs = downwind_deviations_tseries_singfly_all(fly_n, :, 2);
%     xy_speeds = xydists_angles_tseries_all(fly_n, :, 2)./frame_time;
%     [grouped_vals_vec2] = bin_vec2_by_vec1(curr_devs, xy_speeds, angle_bins);   %xyspeeds grouped by corress heading angle for a given fly
%     mean_spd_vecs_all(fly_n, :) = mean(grouped_vals_vec2, 1, 'omitnan');
% end
% mean_spds = mean(mean_spd_vecs_all, 1, 'omitnan');
% se_spds = std(mean_spd_vecs_all, [], 1, 'omitnan')./sqrt(size(mean_spd_vecs_all, 1));
% shadedErrorBar((angle_bins(2:end) - 10), mean_spds, se_spds, {'Color', paired_color}, 1);
% hold on
% 
% %logging grouped speeds for statistical testing
% speeds_ang = pad_n_concatenate(reshape(mean_spd_vecs_all(:, 1:3), [], 1), reshape(mean_spd_vecs_all(:, 7:9), [], 1), 2, nan);
% all_speeds_ang = pad_n_concatenate(all_speeds_ang, speeds_ang, 2, nan);
% 
% 
% %unpaired odor, 25s
% mean_spd_vecs_all = zeros(size(downwind_deviations_tseries_singfly_all, 1), (length(angle_bins)) - 1) + nan;
% for fly_n = 1:size(downwind_deviations_tseries_singfly_all, 1)
%     curr_devs = downwind_deviations_tseries_singfly_all(fly_n, :, 4);
%     xy_speeds = xydists_angles_tseries_all(fly_n, :, 4)./frame_time;
%     [grouped_vals_vec2] = bin_vec2_by_vec1(curr_devs, xy_speeds, angle_bins);   %xyspeeds grouped by corress heading angle for a given fly
%     mean_spd_vecs_all(fly_n, :) = mean(grouped_vals_vec2, 1, 'omitnan');
% end
% mean_spds = mean(mean_spd_vecs_all, 1, 'omitnan');
% se_spds = std(mean_spd_vecs_all, [], 1, 'omitnan')./sqrt(size(mean_spd_vecs_all, 1));
% shadedErrorBar((angle_bins(2:end) - 10), mean_spds, se_spds, {'Color', unpaired_color}, 1);
% hold off
% 
% 
% 
% %logging grouped speeds for statistical testing
% speeds_ang = pad_n_concatenate(reshape(mean_spd_vecs_all(:, 1:3), [], 1), reshape(mean_spd_vecs_all(:, 7:9), [], 1), 2, nan);
% all_speeds_ang = pad_n_concatenate(all_speeds_ang, speeds_ang, 2, nan);
% 
% title('25s gap')
% xlabel('abs. deviation from upwind (degrees)')
% ylabel('xy speed (mm/s)')
% ax_vals = axis;
% ax_vals(1) = 0;
% ax_vals(2) = 180;
% ax_vals(3) = 0;
% ax_vals(4) = 15;
% axis(ax_vals);
% 
% fig_wrapup(36, [], [25, 30], .6);
% 
% 
% %Plotting orientation v/s angular velocity for all time-points in analysis window
% figure(37)
% set(gcf, 'Name', 'downwind deviations v/s angular velocities, 0s')
% %grouping xyspeeds according to their corressponding orientation angles, for each fly
% all_angvels_ang = [];
% %paired odor, 0s
% mean_angvel_vecs_all = zeros(size(downwind_deviations_tseries_all, 1), (length(angle_bins)) - 1) + nan;
% for fly_n = 1:size(downwind_deviations_tseries_all, 1)
%     curr_devs = downwind_deviations_tseries_all(fly_n, :, 1);
%     curr_ang_vels = [nan, abs(movmean(diff(curr_devs), round(0.2/frame_time))./frame_time)];    %frame-by-frame change in orientation followed by boxcar filtering, followed by taking abs val
%     [grouped_angvels_vec2] = bin_vec2_by_vec1(curr_devs, curr_ang_vels, angle_bins);   %xyspeeds grouped by corress heading angle for a given fly
%     mean_angvel_vecs_all(fly_n, :) = mean(grouped_angvels_vec2, 1, 'omitnan');
% end
% 
% 
% mean_angvels = mean(mean_angvel_vecs_all, 1, 'omitnan');                                          %average across flies
% se_angvels = std(mean_angvel_vecs_all, [], 1, 'omitnan')./sqrt(size(mean_angvel_vecs_all, 1));       %SE across flies
% shadedErrorBar((angle_bins(2:end) - 10), mean_angvels, se_angvels, {'Color', paired_color}, 1);   %subtracting 10 from bins to set to middle rather than upper edge of each bin
% hold on
% 
% %logging grouped speeds for statistical testing
% angvels_ang = pad_n_concatenate(reshape(mean_angvel_vecs_all(:, 1:3), [], 1), reshape(mean_angvel_vecs_all(:, 7:9), [], 1), 2, nan);
% all_angvels_ang = pad_n_concatenate(all_angvels_ang, angvels_ang, 2, nan);
% 
% 
% %unpaired odor, 0s
% mean_angvel_vecs_all = zeros(size(downwind_deviations_tseries_all, 1), (length(angle_bins)) - 1) + nan;
% for fly_n = 1:size(downwind_deviations_tseries_all, 1)
%     curr_devs = downwind_deviations_tseries_all(fly_n, :, 3);
%     curr_ang_vels = [nan, abs(movmean(diff(curr_devs), round(0.2/frame_time))./frame_time)];    %frame-by-frame change in orientation followed by boxcar filtering, followed by taking abs val
%     [grouped_angvels_vec2] = bin_vec2_by_vec1(curr_devs, curr_ang_vels, angle_bins);   %xyspeeds grouped by corress heading angle for a given fly
%     mean_angvel_vecs_all(fly_n, :) = mean(grouped_angvels_vec2, 1, 'omitnan');
% end
% mean_angvels = mean(mean_angvel_vecs_all, 1, 'omitnan');                                          %average across flies
% se_angvels = std(mean_angvel_vecs_all, [], 1, 'omitnan')./sqrt(size(mean_angvel_vecs_all, 1));       %SE across flies
% shadedErrorBar((angle_bins(2:end) - 10), mean_angvels, se_angvels, {'Color', unpaired_color}, 1);   %subtracting 10 from bins to set to middle rather than upper edge of each bin
% hold off
% 
% %logging grouped speeds for statistical testing
% angvels_ang = pad_n_concatenate(reshape(mean_angvel_vecs_all(:, 1:3), [], 1), reshape(mean_angvel_vecs_all(:, 7:9), [], 1), 2, nan);
% all_angvels_ang = pad_n_concatenate(all_angvels_ang, angvels_ang, 2, nan);
% 
% 
% title('0s gap')
% xlabel('abs. deviation from upwind (degrees)')
% ylabel('angular vel. (deg/s)')
% ax_vals = axis;
% ax_vals(1) = 0;
% ax_vals(2) = 180;
% ax_vals(3) = -45;
% ax_vals(4) = 45;
% axis(ax_vals);
% 
% fig_wrapup(37, [], [25, 30], .6);
% 
% figure(38)
% set(gcf, 'Name', 'downwind deviations v/s angular velocities, 25s')
% %grouping xyspeeds according to their corressponding orientation angles, for each fly
% all_angvels_ang = [];
% %paired odor, 25s
% mean_angvel_vecs_all = zeros(size(downwind_deviations_tseries_all, 1), (length(angle_bins)) - 1) + nan;
% for fly_n = 1:size(downwind_deviations_tseries_all, 1)
%     curr_devs = downwind_deviations_tseries_all(fly_n, :, 2);
%     curr_ang_vels = [nan, abs(movmean(diff(curr_devs), round(0.2/frame_time))./frame_time)];    %frame-by-frame change in orientation followed by boxcar filtering, followed by taking abs val
%     [grouped_angvels_vec2] = bin_vec2_by_vec1(curr_devs, curr_ang_vels, angle_bins);   %xyspeeds grouped by corress heading angle for a given fly
%     mean_angvel_vecs_all(fly_n, :) = mean(grouped_angvels_vec2, 1, 'omitnan');
% end
% 
% 
% mean_angvels = mean(mean_angvel_vecs_all, 1, 'omitnan');                                          %average across flies
% se_angvels = std(mean_angvel_vecs_all, [], 1, 'omitnan')./sqrt(size(mean_angvel_vecs_all, 1));       %SE across flies
% shadedErrorBar((angle_bins(2:end) - 10), mean_angvels, se_angvels, {'Color', paired_color}, 1);   %subtracting 10 from bins to set to middle rather than upper edge of each bin
% hold on
% 
% %logging grouped speeds for statistical testing
% angvels_ang = pad_n_concatenate(reshape(mean_angvel_vecs_all(:, 1:3), [], 1), reshape(mean_angvel_vecs_all(:, 7:9), [], 1), 2, nan);
% all_angvels_ang = pad_n_concatenate(all_angvels_ang, angvels_ang, 2, nan);
% 
% 
% %unpaired odor, 25s
% mean_angvel_vecs_all = zeros(size(downwind_deviations_tseries_all, 1), (length(angle_bins)) - 1) + nan;
% for fly_n = 1:size(downwind_deviations_tseries_all, 1)
%     curr_devs = downwind_deviations_tseries_all(fly_n, :, 4);
%     curr_ang_vels = [nan, abs(movmean(diff(curr_devs), round(0.2/frame_time))./frame_time)];    %frame-by-frame change in orientation followed by boxcar filtering, followed by taking abs val
%     [grouped_angvels_vec2] = bin_vec2_by_vec1(curr_devs, curr_ang_vels, angle_bins);   %xyspeeds grouped by corress heading angle for a given fly
%     mean_angvel_vecs_all(fly_n, :) = mean(grouped_angvels_vec2, 1, 'omitnan');
% end
% mean_angvels = mean(mean_angvel_vecs_all, 1, 'omitnan');                                          %average across flies
% se_angvels = std(mean_angvel_vecs_all, [], 1, 'omitnan')./sqrt(size(mean_angvel_vecs_all, 1));       %SE across flies
% shadedErrorBar((angle_bins(2:end) - 10), mean_angvels, se_angvels, {'Color', unpaired_color}, 1);   %subtracting 10 from bins to set to middle rather than upper edge of each bin
% hold off
% 
% %logging grouped speeds for statistical testing
% angvels_ang = pad_n_concatenate(reshape(mean_angvel_vecs_all(:, 1:3), [], 1), reshape(mean_angvel_vecs_all(:, 7:9), [], 1), 2, nan);
% all_angvels_ang = pad_n_concatenate(all_angvels_ang, angvels_ang, 2, nan);
% 
% 
% title('25s gap')
% xlabel('abs. deviation from upwind (degrees)')
% ylabel('angular vel. (deg/s)')
% ax_vals = axis;
% ax_vals(1) = 0;
% ax_vals(2) = 180;
% ax_vals(3) = -45;
% ax_vals(4) = 45;
% axis(ax_vals);
% 
% fig_wrapup(38, [], [25, 30], .6);
% 
% %Doing statistics on xy speeds, separated by orientation angle
% figure(39)
% %0s gap
% set(gcf, 'Name', 'xyspeeds grouped by angles, statsitics, 0s')
% score_vecs_all_final = all_speeds_ang(:, 1:4);        
% markercolor = [paired_color; paired_color; unpaired_color; unpaired_color];
% xlabels = [{'prd, 0-60'}, {'prd, 120-180 deg'}, {'unprd, 0-60'}, {'unprd, 120-180 deg'}];
% [fig_h, r_vecs_saved] = scattered_dot_plot_ttest(score_vecs_all_final, 39, 2.5, 4, 4, markercolor, 1, [], [], xlabels, 2, [0, 0, 0], 2, 0.05, 0, 1, 'force_mean');
% title('xyspeeds grouped by angles, 0s')
% ax_vals = axis;
% ax_vals(4) = 30;
% axis(ax_vals);
% ylabel('speed (mm/s)');
% fig_wrapup(fig_h, [], [25, 30], .6);
% 
% figure(40)
% %25s gap
% set(gcf, 'Name', 'xyspeeds grouped by angles, statsitics, 25s')
% score_vecs_all_final = all_speeds_ang(:, 5:8);        
% markercolor = [paired_color; paired_color; unpaired_color; unpaired_color];
% xlabels = [{'prd, 0-60'}, {'prd, 120-180 deg'}, {'unprd, 0-60'}, {'unprd, 120-180 deg'}];
% [fig_h, r_vecs_saved] = scattered_dot_plot_ttest(score_vecs_all_final, 40, 2.5, 4, 4, markercolor, 1, [], [], xlabels, 2, [0, 0, 0], 2, 0.05, 0, 1, 'force_mean');
% title('xyspeeds grouped by angles, 25s')
% ax_vals = axis;
% ax_vals(4) = 30;
% axis(ax_vals);
% ylabel('speed (mm/s)');
% fig_wrapup(fig_h, [], [25, 30], .6);
% 
% % %Doing statistics on xy speeds, separated by orientation angle
% % figure(41)
% % %0s gap
% % set(gcf, 'Name', 'ang vels grouped by angles, statsitics, 0s')
% % score_vecs_all_final = all_angvels_ang(:, 1:4);        
% % markercolor = [paired_color; paired_color; unpaired_color; unpaired_color];
% % xlabels = [{'prd, 0-60'}, {'prd, 120-180 deg'}, {'unprd, 0-60'}, {'unprd, 120-180 deg'}];
% % [fig_h, r_vecs_saved] = scattered_dot_plot_ttest(score_vecs_all_final, 41, 2.5, 4, 4, markercolor, 1, [], [], xlabels, 1, [0, 0, 0], 2, 0.05, 0);
% % title('angvels grouped by angles, 0s')
% % ax_vals = axis;
% % ax_vals(4) = 30;
% % axis(ax_vals);
% % ylabel('ang. vel. (deg./s)');
% % fig_wrapup(fig_h, [], [25, 30], .6);
% % 
% % figure(42)
% % %25s gap
% % set(gcf, 'Name', 'ang vels grouped by angles, statsitics, 25s')
% % score_vecs_all_final = all_angvels_ang(:, 5:8);        
% % markercolor = [paired_color; paired_color; unpaired_color; unpaired_color];
% % xlabels = [{'prd, 0-60'}, {'prd, 120-180 deg'}, {'unprd, 0-60'}, {'unprd, 120-180 deg'}];
% % [fig_h, r_vecs_saved] = scattered_dot_plot_ttest(score_vecs_all_final, 42, 2.5, 4, 4, markercolor, 1, [], [], xlabels, 1, [0, 0, 0], 2, 0.05, 0);
% % title('angvels grouped by angles, 25s')
% % ax_vals = axis;
% % ax_vals(4) = 30;
% % axis(ax_vals);
% % ylabel('ang. vel. (deg./s)');
% % fig_wrapup(fig_h, [], [25, 30], .6);
% % 
% 

keyboard


%plotting trajectory video
if plot_video == 1
    more_frames = 1;
   
    for dset_type = 1:4

        win_width = round( (t_window_orig(2) - t_window_orig(1))./frame_time);
        if dset_type == 1 | dset_type == 3  %0s gap case
            f0 = round( (t_window_orig(1) + pulse_times(1))./frame_time );
        elseif dset_type == 2 | dset_type == 4  %25s gap case
            f0 = round( (t_window_orig(1) + pulse_times(2))./frame_time );
        else
        end
        f1 = f0 - round(win_width./4);
        f2 = f0 + win_width;

        %computing time-stamp vector
        neg_vec = flipud((0:-frame_time:-((win_width./4).*frame_time))');    
        pos_vec = (0:frame_time:win_width.*frame_time)';                     
        t_vec = round([neg_vec; pos_vec(2:end)].*1000);     %in ms

        traj_mat = traj_mat_exps_all(:, f1:f2, :, dset_type);

        if align_vids_vert == 1
            %subtracting away t0 angular offset to align all flies to the upward vertical on frame1
            [thetaxy, rho] = cart2pol(traj_mat(:, :, 1), traj_mat(:, :, 2));
            thetaxy1 = thetaxy(:, round(win_width./4));
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
        
       %assigning a different variable name to each set of adjusted co-ords
       eval(['incl_flies', num2str(dset_type), '= squeeze(traj_mat_exps_all(:, 1, 1, dset_type));']);
       eval(['traj_mat', num2str(dset_type), '= traj_mat;']);
       
    end
        
    set_names = [{'no gap, unpaired'}, {'25s gap, unpaired'}, {'no gap, paired'}, {'25s gap, paired'}];
    
    figure(99);
    v_all = VideoWriter([vid_path, 'MBON_act'], 'MPEG-4');
    open(v_all)
    red_vec = (0:50)'./50;
    color_vec = [red_vec, zeros(51, 2)];
    
    for frame_n = 1:(size(traj_mat, 2))

        for dset_n = 1:4
            if dset_n == 2
                subplot_tight(2, 2, 3, 0.05)
            elseif dset_n == 3
                subplot_tight(2, 2, 2, 0.05)
            else
                subplot_tight(2, 2, dset_n, 0.05)
            end
                
                
            set_name = set_names{dset_n};
            eval(['traj_mat = traj_mat', num2str(dset_n), ';']);
            
    
            for fly_n = 1:size(traj_mat, 1)
                curr_traj = traj_mat(fly_n, 1:end, :);
                
                eval(['incl_flies = incl_flies', num2str(dset_n), ';']);
                if incl_flies(fly_n, 1) == 0       %case where fly was excluded for being close to the edge/center
                    %use_color = [0.95, 0.95, 0.95];
                    continue
%                 elseif incl_flies(fly_n, 1) == 1   %case where fly was included
%                     use_color = [0, 0, 0];
                else
                end
                if isnan(curr_traj(1, frame_n, 1)) == 0
                    %computing curr_dist
                    curr_dist = mean(sqrt(sum(curr_traj(:, frame_n, 1:2).^2, 3)), 2, 'omitnan');
                    use_color = color_vec(round(curr_dist), :);
                    if curr_dist > 51
                        continue
                    else
                    end

                    [pt2] = get_orientation_pt([curr_traj(1, frame_n, 1), curr_traj(1, frame_n, 2)], curr_traj(1, frame_n, 3), .5); %using tracked orientation to specify arrow direction
                    arr_pts = [curr_traj(1, (frame_n), 1), curr_traj(1, (frame_n), 2), pt2(1), pt2(2)];
                    plot_arrow(arr_pts(1), arr_pts(2), arr_pts(3), arr_pts(4), 'color', use_color, 'head_width', 2, 'head_height', 3, 'rel_sizing', 0, 'EdgeColor', [1, 1, 1], 'FaceColor', use_color);
                    if dset_n == 4
                        t = text(35, -45, [num2str(t_vec(frame_n)), ' ms']);
                        t.FontSize = 12;
                        if t_vec(frame_n) >=0
                            t.Color = 'red';
                        else
                        end
                    else
                    end
                    if frame_n == 1
                        title(set_names{dset_n});
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
            ax = gca;
            
        end
        hold off
        %writing fig frame to video file
        set(gcf, 'Position', [305, 136, 840, 630]);
        frame = getframe(gcf);
        cla(ax);
        writeVideo(v_all,frame);
        
    end
    close(v_all);
    
    
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

function [traj_samps, upwind_dists, traj_mat_ext, upwind_dists_tseries, traj_mat_exp, radial_pos_tseries] = compute_center_dists(traj_mat, frame_time, r_cutoff, normalize_cent_dists, equilib_time, pulse_width, win_width, initial_delay, use_occupancy)
    
    %saving all trajectories, but logging which flies were not discarded as edge flies
    traj_mat_exp = zeros(size(traj_mat, 1), (size(traj_mat, 2) + 1), size(traj_mat, 3)) + 1;
    traj_mat_exp(:, 2:end, :) = traj_mat;
        
    traj_mat_ext = traj_mat;    
    traj_mat = traj_mat(:, :, 1:2);
    
    frame0 = round(initial_delay./frame_time); 
    loc_0 = squeeze(traj_mat(:, frame0, :));     %location at beginning of t_window, all flies
    edge_fliesi = find(sqrt(sum(loc_0.^2, 2)) < r_cutoff(1) | sqrt(sum(loc_0.^2, 2)) > r_cutoff(2));
    traj_mat(edge_fliesi, :, :) = [];               %excluding flies too close to edge or center
    traj_mat_ext(edge_fliesi, :, :) = [];
    traj_mat_exp(edge_fliesi, 1, :) = 0;
    
    upwind_dists_tseries = [];
    radial_pos_tseries = [];
    fnum = 1;
    
    win_wing = round(((win_width - pulse_width)./2)./frame_time);    %the trajectory sampling time window on either side of LED pulse delivery time
    for frame_n = (round(initial_delay./frame_time) - win_wing):1:(round((initial_delay + pulse_width)./frame_time) + win_wing)
        if fnum == 1
            try
                zero_pts = squeeze(traj_mat(:, frame_n, :));  %t0 is always the time point of the transition
            catch
                keyboard
            end
            
            zero_dists = sqrt(sum(zero_pts.^2, 2));         %distance from 0 at t0
            if normalize_cent_dists == 1
                zero_dists = (zero_dists./50).^2;
            else
            end
        else
        end
        
        curr_pts = squeeze(traj_mat(:, frame_n, :));
        curr_dists_abs = sqrt(sum(curr_pts.^2, 2));         %distance from 0 in current frame
       
        %normalizing distances from center by area if manually specified
        if normalize_cent_dists == 1
            curr_dists_abs = (curr_dists_abs./50).^2;
        else
        end
        curr_dists = curr_dists_abs - zero_dists;           %computing delta dist relative to t0
        upwind_dists_tseries = [upwind_dists_tseries, curr_dists];
        radial_pos_tseries = [radial_pos_tseries, curr_dists_abs];
        fnum = fnum + 1;
        
    end
    
    if use_occupancy == 0
        upwind_dists = upwind_dists_tseries(:, (round(pulse_width./frame_time) + win_wing) );    %upwind distance at time point at end of LED pulse, but in sampled traces
    elseif use_occupancy == 1
        upwind_dists = mean(upwind_dists_tseries(:, win_wing:(round(pulse_width./frame_time) + win_wing)), 2, 'omitnan');      %mean upwind distance over entire LED pulse duration
    else
    end
   
   
    traj_samps = traj_mat(:, (round(initial_delay./frame_time) - win_wing):(round((initial_delay + pulse_width)./frame_time) + win_wing), :);
    
    %identifying and removing flies that don't move throughout the analysis time window
%     [xy_vels_tseries, xy_vels_bin, ststp_tseries] = get_xydists_ststp_probs(traj_mat, frame_time, t_window, [2, 0.5]);
%     tot_moving_frames = sum(xy_vels_bin, 2, 'omitnan');
%     non_movers = find(tot_moving_frames == 0);
%     traj_mat(non_movers, :, :) = [];               %excluding flies that aren't moving
%     traj_mat_ext(non_movers, :, :) = [];
%     traj_mat_exp(non_movers, 1, :) = 0;
%     traj_samps(non_movers, :, :) = [];
    
end


function [mean_downwind_deviations, downwind_deviations_tseries, edge_flies, xy_dists_tseries] = compute_radial_orientations(traj_mat, frame_time, t_window, ang_cutoff, edge_cutoff, analysis_offset)
    
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
    ori_mat(del) = 360 - abs(ori_mat(del)) + pi./2; 
    
    %computing deviation from upwind orientation (0 = perfectly upwind orientation)
    radial_orientations_mat = ori_mat - theta;
    
    %correcting for abs(angles) > 180
    del = find(radial_orientations_mat > 180);
    radial_orientations_mat(del) = radial_orientations_mat(del) - 360;
    del = find(radial_orientations_mat < -180);
    radial_orientations_mat(del) = radial_orientations_mat(del) + 360;
    
    radial_orientations_mat = abs(radial_orientations_mat);
    
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
    t_window_curr = t_window - analysis_offset;
    fr1 = max([round(t_window_curr(1)./frame_time), 1]);
    fr2 = round(t_window_curr(2)./frame_time);
    
    downwind_deviations_tseries = radial_orientations_mat(:, fr1:fr2);
    
    %identifying flies <2mm from upwind edge to highlight during plotting
    edge_flies = zeros(size(traj_mat, 1), 1);
    dists = mean(sqrt(sum(traj_mat(:, fr1:fr2, 1:2).^2, 3)), 2, 'omitnan');
    edge_flies(dists>edge_cutoff) = 1;
    
    %computing mean downind deviation
    mean_downwind_deviations = mean(abs(downwind_deviations_tseries), 2, 'omitnan');
    
    %computing xy dist travelled in each frame (ie. xy speed)
    traj_mat_t_win = traj_mat(:, (fr1 - 1):fr2, :);
    for frame_n = 2:size(traj_mat_t_win, 2)
        if frame_n > 2
            pos1 = pos0;
        else
            pos1 = traj_mat_t_win(:, frame_n, 1:2);
        end
        pos0 = traj_mat_t_win(:, (frame_n - 1), 1:2);
        
        %computing distances
        xy_dists_tseries(:, (frame_n - 1)) = sqrt((pos1(:, :, 1) - pos0(:, :, 1)).^2 + (pos1(:, :, 2) - pos0(:, :, 2)).^2);
    end
    
end

function [mean_abs_dists] = compute_abs_dists(traj_mat, frame_time, t_window, analysis_offset)
    %sampling tseries in t_window
    t_window_curr = t_window - analysis_offset;
    fr1 = round(t_window_curr(1)./frame_time);
    fr2 = round(t_window_curr(2)./frame_time);
    mean_abs_dists = mean(sqrt(sum(traj_mat(:, fr1:fr2, 1:2).^2, 3)), 2, 'omitnan');
end

function [pt2] = get_orientation_pt(pt1, angle, dist)
    %This function uses a point in xy space, and an orientation angle in
    %radians to specify a second point in xy space that is at dist and angle
    %relative to the first point.
    [x0, y0] = pol2cart(angle, dist);     %pt2 at specified angle and distance relative to the origin
    pt2 = [x0, y0] + pt1;

end

function [xy_vels_tseries, xy_vels_bin, ststp_tseries] = get_xydists_ststp_probs(traj_mat, frame_time, t_window, vel_cutoffs, analysis_offset)

    t_window_curr = t_window - analysis_offset;
    t_window = round(t_window_curr./frame_time);
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

