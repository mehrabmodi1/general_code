function [sig_cell_mat, all_bad_trs] = cell_data_quality_control(dff_data_mat, stim_mat, stim_mat_simple, column_heads, sig_cell_mat, sig_cell_mat_key, suppress_plots, frame_time)
%Syntax: [sig_cell_mat, all_bad_trs] = cell_data_quality_control(dff_data_mat, stim_mat, stim_mat_simple, column_heads, sig_cell_mat, manual_inspec)
%This function plots a number of useful visualisations of reliability of
%cell responses. It forces to 0, cells in sig_cell_mat that were found
%significant, but show > 25% variability in more than half the trials.
%all_bad_trs is a list of trials with > 25% variability in the
%cell-averaged responses.

olf1_od_coln = find_stim_mat_simple_col('odor_n', column_heads);
olf2_od_coln = find_stim_mat_simple_col('odour_olf2', column_heads);
olf1_dur_coln = find_stim_mat_simple_col('duration', column_heads);
olf2_dur_coln = find_stim_mat_simple_col('duration_olf2', column_heads);


global color_vec
global greymap

all_bad_trs = [];
for stim_gp_n = 1:size(sig_cell_mat_key, 1)
    curr_stim_gp = sig_cell_mat_key(stim_gp_n, :);      %stimulus parameters for current trial group
    odor_ni = stim_gp_n;
    curr_trs = find(stim_mat_simple(:, olf1_od_coln) == curr_stim_gp(1) & stim_mat_simple(:, olf1_dur_coln) == curr_stim_gp(2) &...
                                                    stim_mat_simple(:, olf2_od_coln) == curr_stim_gp(3) & stim_mat_simple(:, olf2_dur_coln) == curr_stim_gp(4));
    
    od_pulse_frames = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
    
    olfn_used = find(curr_stim_gp(1, [2, 4]) > 0);
    od_pulse_frames = od_pulse_frames{olfn_used(1)};
    
    ave_mat = mean(dff_data_mat(:, :, curr_trs), 3, 'omitnan');
    try
        curr_sig_cells = sig_cell_mat(:, stim_gp_n);
    catch
        keyboard
    end
    
    try
        curr_color = color_vec(odor_ni, :);
    catch
        curr_color = [0.6, 0.6, 0.6];
    end
            
    
    %making plots for manual inspection only - to check if sig cells have been found or left out
    if suppress_plots == 0
        curr_sig_cells1 = find(curr_sig_cells == 1);
        curr_sig_cells0 = find(curr_sig_cells == 0);
        max_val = 5;
        figure(1)
        imagesc(ave_mat(:, curr_sig_cells1)', [0, max_val])
        colormap(greymap)
        ylabel('significant cell number')
        set_xlabels_time(1, frame_time, 10);
        fig_wrapup(1, [])
        try
            
            add_stim_bar(1, od_pulse_frames, curr_color);
        catch
            keyboard
        end
        
        figure(2)
        imagesc(ave_mat(:, curr_sig_cells0)', [0, max_val])
        colormap(greymap)
        ylabel('non-significant cell number')
        set_xlabels_time(2, frame_time, 10);
        fig_wrapup(2, [])
        add_stim_bar(2, od_pulse_frames, curr_color);

    else
    end
    curr_sig_cells = find(curr_sig_cells == 1);    

    %checking if responses in first-last half of trials are similar to each other and to overall mean.
    q_n = floor(length(curr_trs)./2);

    if q_n > 0
        ave_mat_early = mean(dff_data_mat(:, curr_sig_cells, curr_trs(1:q_n)), 3, 'omitnan');
        ave_mat_late = mean(dff_data_mat(:, curr_sig_cells, curr_trs(((q_n) + 1): end)), 3, 'omitnan');

        curr_sig_cells = sig_cell_mat(:, stim_gp_n);
        curr_sig_cells = find(curr_sig_cells == 1);

        if suppress_plots == 0
            fig_h = figure(3);
            imagesc(ave_mat_early', [0, max_val])
            colormap(greymap)
            ylabel('significant cell number')
            set_xlabels_time(3, frame_time, 10);
            fig_wrapup(3, [])
            add_stim_bar(3, od_pulse_frames, curr_color);
            set(fig_h, 'Name','ave resp traces for first half of trials')

            fig_h = figure(4);
            imagesc(ave_mat_late', [0, max_val])
            colormap(greymap)
            ylabel('significant cell number')
            set_xlabels_time(4, frame_time, 10);
            fig_wrapup(4, [])
            add_stim_bar(4, od_pulse_frames, curr_color);
            set(fig_h, 'Name','ave resp traces for last half of trials')


        else
        end



    else
        disp('too few reps to analyse by quarters of reps... skipping.')
    end

    %computing difference between matrices
    %using ave_mat as the reference matrix and calculating
    %corrcoefs as a measure of similarity
    corr_mat = zeros(length(curr_sig_cells), length(curr_trs));
    for rep_n = 1:length(curr_trs)
        curr_tr = curr_trs(rep_n);
        curr_traces = dff_data_mat(:, curr_sig_cells, curr_tr);     %single trial response matrix of size n_frames by n_sig_cells
        for cell_n = 1:length(curr_sig_cells)
            ref_trace = ave_mat(:, curr_sig_cells(cell_n) );
            curr_trace = curr_traces(:, cell_n);
            corr = corrcoef(ref_trace, curr_trace);
            corr_mat(cell_n, rep_n) = corr(1, 2);
        end

    end

    if suppress_plots == 0
        figure(5)
        imagesc(corr_mat, [0, 1])
        xlabel('curr trial n')
        ylabel('sig cell n')
        title('corr w mean resp vector')

        figure(6)
        plot(mean(corr_mat))
        ylabel('cell-averaged corr w mean resp vector')
        xlabel('curr trial n')
        ax_vals = axis;
        ax_vals(1, [3,4]) = [0, 1];
        axis(ax_vals);
    else
    end

    %checking if template-corr averaged across all cells is flat
    %across trials and logging a trial as bad if if variation is > 25%
    mean_corrs = mean(corr_mat, 1, 'omitnan');
    bad_trs = find(mean_corrs < 0.75.*max(mean_corrs));
    mean_corrs(:, bad_trs) = mean_corrs(:, bad_trs) + nan;
    bad_trs = curr_trs(bad_trs);
    all_bad_trs = [all_bad_trs; bad_trs];

    %checking for big swings in template-corr for each cell and
    %throwing cells away if swing is > 25% for more than half the
    %trials
    n_good_reps = size(mean_corrs, 2) - sum(isnan(mean_corrs));
    bad_cells = [];
    for sig_cell_n = 1:size(corr_mat, 1)
        curr_bad_trs = find(corr_mat(sig_cell_n, :) < 0.75.*max(corr_mat(sig_cell_n, :)));

        if length(curr_bad_trs) > n_good_reps./2
            sig_cell_ni = curr_sig_cells(sig_cell_n);
            bad_cells = [bad_cells; sig_cell_ni];

        else
        end

    end

    %removing bad cells from sig_cell_mat
    sig_cell_mat(bad_cells, odor_ni) = 0;



     if suppress_plots == 0
         %checking again if responses in first-last half of trials are similar to each other and to overall mean.
         q_n = floor(length(curr_trs)./2);

        if q_n > 0
            curr_sig_cells = sig_cell_mat(:, odor_ni);
            curr_sig_cells = find(curr_sig_cells == 1);

            ave_mat_early = mean(dff_data_mat(:, curr_sig_cells, curr_trs(1:q_n)), 3, 'omitnan');
            ave_mat_late = mean(dff_data_mat(:, curr_sig_cells, curr_trs(((q_n) + 1): end)), 3, 'omitnan');

            fig_h = figure(7);
            imagesc(ave_mat_early', [0, max_val])
            colormap(greymap)
            ylabel('significant cell number')
            set_xlabels_time(7, frame_time, 10);
            fig_wrapup(7, [])
            add_stim_bar(7, od_pulse_frames, curr_color);
            set(fig_h, 'Name','ave resp traces for first half of trials')

            fig_h = figure(8);
            imagesc(ave_mat_late', [0, max_val])
            colormap(greymap)
            ylabel('significant cell number')
            set_xlabels_time(8, frame_time, 10);
            fig_wrapup(8, [])
            add_stim_bar(8, od_pulse_frames, curr_color);
            set(fig_h, 'Name','ave resp traces for last half of trials')


        else
        end
        
        
        %plotting resps for each trial
        figure(9)
        curr_traces = dff_data_mat(:, :, curr_trs);
        curr_resps = mean(curr_traces(od_pulse_frames(1):od_pulse_frames(2), :, :), 1, 'omitnan');
        curr_sig_cells = sig_cell_mat(:, stim_gp_n);
        curr_sig_cells = find(curr_sig_cells == 1);
        imagesc(squeeze(curr_resps(1, curr_sig_cells, :))');
        colormap(greymap)
        xlabel('cell number')
        ylabel('trial number')
        keyboard
    else
    end

    close all
    
end