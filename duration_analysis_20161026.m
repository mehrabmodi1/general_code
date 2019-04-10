clear all
close all

direc_lists_mat = [{'D:\Data\CSHL\dataset_list_sustained_MB418B_20160414.txt'};... %KC A'B'
                    {'D:\Data\CSHL\dataset_list_sustained_MB185B_20160426.txt'};... %KC AB
                    %{'D:\Data\CSHL\dataset_list_stim_dur_20160316.txt'} ... KCs all
                    {'D:\Data\CSHL\dataset_list_sustained_MB131B_20160515.txt'}... %KC G
                    ]; 
                
dump_direcs = [{'C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Analysis\data_dump_20160624\ApBp\'};...
                {'C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Analysis\data_dump_20160624\AB\'};...
                {'C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Analysis\data_dump_20160624\G\'}];
            
            
n_direc_lists = size(direc_lists_mat, 1);
                
                
color_vec = load('C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Code\std_color_vec.txt');
a = colormap('bone');
greymap = flipud(a);

suppress_plots = 1;       %1 - doesn't plot stuff, 0 - plots stuff

saved_cell_data = cell(1, 45);

%figure property initialisation variables
plot_height = 200;
plot_width = 200;
axis_font_size = 8;
tick_length = [0.02, 0.05];
line_widths = 2;
marker_sizes = 1;
marker_sizes_f = 4;
box_l_width = 0.5;


saved_resp_vecs_all = [];
saved_resp_vecs_1s_all = [];


%loop to go through all directory lists
for direc_list_n = 1:n_direc_lists

    list_direc = direc_lists_mat{direc_list_n, 1};
    fid = fopen(list_direc);
    direc_counter = 0;
    saved_resp_vecs = [];
    saved_resp_vecs_1s = [];
    saved_resp_vec_ratio = [];
    all_sig_cells_mat = [];
    all_sig_cells_1s_mat = [];
    
    %loop to go through all experiment datasets listed in list file
    while 1
        direc_counter = direc_counter + 1;
        direc = fgetl(fid);

        if ischar(direc) ~= 1
            break
        else
        end

        %skipping a particular OK107 dataset where n_frames was < stim time
        del = strfind(direc, '20160310');
        if isempty(del) == 0
            continue
        else
        end


        direc = [direc, '\'];

        %replacing C: with D:
        a = strfind(direc, 'C:');
        direc(a) = 'D';

        %skipping if this dataset has already been analysed
    %     if exist([direc 'cell_classifications.mat']) == 2
    %         disp('dataset already analysed; skipping...')
    %         %continue
    %     else
    %     end


        %loading extracted raw fluorescence data matrices written by
        %raw_dff_extractor
        dataset = load([direc 'expt.mat']);
        dataset = dataset.data;
        raw_data_mat = load([direc 'expt_raw_traces.mat']);
        raw_data_mat = raw_data_mat.raw_data_mat;


        %calculating dF/F traces and creating the sparse, 4-D, nan-filled
        %dff_data_mat 
        [dff_data_mat, stim_mat, prot_switch_trials] = cal_dff_traces_20160317(raw_data_mat, dataset, list_direc);
        
        dump_direc = dump_direcs{direc_list_n};
        save([dump_direc 'dFF_data_fly' int2str(direc_counter) '.mat'], 'dff_data_mat');
        save([dump_direc 'stim_info_fly' int2str(direc_counter) '.mat'], 'stim_mat');
        
        clear raw_data_mat
        odor_list = unique(stim_mat(:, 1) );
        del = find(odor_list == 0);
        odor_list(del) = [];
        n_odors = length(odor_list);
        n_frames = size(dff_data_mat, 1);
        n_cells = size(dff_data_mat, 2);
        n_trials = size(dff_data_mat, 3);
        prot_switch_trials = n_trials + 1;                  
        stim_time = dataset(1).stim.stimLatency + 0.54;     %odor stim delay measured with PID is 540 ms
        frame_time = dataset(1).info.framePeriod;
        stim_frame = floor(stim_time./frame_time);
        an_trial_window = nan;
        odor_dur_list = unique(stim_mat(:, 2) );
        n_reps = dataset(1).stim.reps;
        stim_end_frames = [round(odor_dur_list./frame_time) + stim_frame];
        
        
        %identifying sig responses on a single trial basis, and then sig
        %responder cells in any individual block
        [resp_areas, sig_trace_mat, sig_cell_mat] = ...
            cal_sig_responses_20161024(dataset, dff_data_mat, stim_mat, prot_switch_trials, list_direc, an_trial_window);
            
        del = isnan(resp_areas(1, :));
        bad_tr_list = find(del == 1);                   %list of trials thrown away due to movement
        good_tr_list = 1:n_trials;
        good_tr_list(bad_tr_list) = [];
        n_prot_blocks = length(prot_switch_trials);     %number of protocol blocks in this dataset

        %making sure that at least half the trials in each block, for each odor were acquired
        curr_window = n_trials./(length(odor_list).*n_prot_blocks);
        bad_od_blk_counter = 0;
        for odor_n = 1:length(odor_list)
            odor_ni = odor_list(odor_n);
            for odor_dur_n = 1:length(odor_dur_list)
                odor_dur_ni = odor_dur_list(odor_dur_n);
                for block_n = 1:n_prot_blocks
                    od_blk_window_trs = odor_trial_list_builder_20160317(stim_mat, prot_switch_trials, odor_ni, odor_dur_ni, block_n, an_trial_window, 1);

                    n_bad_trs = length(intersect(od_blk_window_trs, bad_tr_list));            %number of bad trials in current block
                    if n_bad_trs > curr_window/2
                        disp('More than half the trials missing for one odor for one block.')
                        bad_od_blk_counter = bad_od_blk_counter + 1;
                    else
                    end

                end
            end
        end
        if bad_od_blk_counter > 1

            %continue
        else
        end
        
        sig_cell_1s_mat = squeeze(sig_cell_mat(:, :, 1));
        sig_cell_mat = squeeze(sig_cell_mat(:, :, 3));
        sig_cell_mat_old = sig_cell_mat;
        
       
        sig_cell_1s_mat_old = sig_cell_1s_mat;
        
        sig_cell_mat = sig_cell_mat(:, odor_list);
        sig_cell_1s_mat = sig_cell_1s_mat(:, odor_list);

        %Keeping track of fraction of responsive cells for each odor
        all_sig_cells_mat = [all_sig_cells_mat; sig_cell_mat];
        all_sig_cells_1s_mat = [all_sig_cells_1s_mat; sig_cell_1s_mat];

        %keeping track of resp cell fractions for each fly
        resp_frac_vec = sum(sig_cell_mat)./size(sig_cell_mat, 1);
        resp_frac_vec_1s = sum(sig_cell_1s_mat)./size(sig_cell_1s_mat, 1);

        saved_resp_vecs = [saved_resp_vecs; resp_frac_vec];
        saved_resp_vecs_1s = [saved_resp_vecs_1s; resp_frac_vec_1s];
        
        saved_resp_vec_ratio = [saved_resp_vec_ratio; (resp_frac_vec - resp_frac_vec_1s)./resp_frac_vec];
        
        
        %looking at ave response traces for responders of various kinds
        sig_mat_diff = sig_cell_mat - sig_cell_1s_mat;
        for cell_n = 1:n_cells
            is_long_responder = sum(sig_mat_diff(cell_n, :));
            if is_long_responder == 0
                continue
            else
            end
            
            sig_odors = find(sig_mat_diff(cell_n, :) == 1);
            if suppress_plots == 0

                for sig_odor_n = 1:is_long_responder
                    sig_odor = sig_odors(sig_odor_n);

                    %list of short stim trials for current cell, current sig_odor
                    short_trs = odor_trial_list_builder_20160317(stim_mat, prot_switch_trials, sig_odor, odor_dur_list(1), block_n, an_trial_window, 1);
                    ave_short_trace = nanmean(squeeze(dff_data_mat(:, cell_n, short_trs, sig_odor)), 2);

                    %list of long stim trials for current cell, current sig_odor
                    long_trs1 = odor_trial_list_builder_20160317(stim_mat, prot_switch_trials, sig_odor, odor_dur_list(2), block_n, an_trial_window, 1);
                    ave_long_trace1 = nanmean(squeeze(dff_data_mat(:, cell_n, long_trs1, sig_odor)), 2);
                    long_trs2 = odor_trial_list_builder_20160317(stim_mat, prot_switch_trials, sig_odor, odor_dur_list(3), block_n, an_trial_window, 1);
                    ave_long_trace2 = nanmean(squeeze(dff_data_mat(:, cell_n, long_trs2, sig_odor)), 2);

                    h_fig = figure(3);
                    plot(ave_short_trace, 'k', 'LineWidth', 2)
                    add_stim_shading(3, [stim_frame, stim_end_frames(1)], 0.5, color_vec(sig_odor, :))
                    set_xlabels_time(3, frame_time, 1)
                    xlabel('time (s)')
                    ylabel('dF/F')
                    set(h_fig, 'units','normalized','position',[.1 .1 .25 .35])
                    
                    h_fig = figure(4);
                    plot(ave_long_trace1, 'k', 'LineWidth', 2)
                    add_stim_shading(4, [stim_frame, stim_end_frames(2)], 0.5, color_vec(sig_odor, :))
                    set_xlabels_time(4, frame_time, 0.6)
                    xlabel('time (s)')
                    ylabel('dF/F')
                    set(h_fig, 'units','normalized','position',[.2 .2 .25 .35])
                    
                    h_fig = figure(5);
                    plot(ave_long_trace2, 'k', 'LineWidth', 2)
                    add_stim_shading(5, [stim_frame, stim_end_frames(3)], 0.5, color_vec(sig_odor, :))
                    set_xlabels_time(5, frame_time, 0.4)
                    xlabel('time (s)')
                    ylabel('dF/F')
                    set(h_fig, 'units','normalized','position',[.3 .3 .25 .35])
                    

                    del = input('Next cell?');
                end
            else
            end
            
        end
       
        
    end
    fclose(fid);

    saved_resp_vecs_all = pad_n_concatenate(saved_resp_vecs_all, saved_resp_vecs, 3, nan);              %saving each fly's resp vec across dataset lists
    saved_resp_vecs_1s_all = pad_n_concatenate(saved_resp_vecs_1s_all, saved_resp_vecs_1s, 3, nan);     %saving each fly's resp vec across dataset lists
    
    %plotting distribution of responsive cells for each odor
    n_resp_vec = sum(all_sig_cells_mat)./size(all_sig_cells_mat, 1);
    n_resp_vec_1s = sum(all_sig_cells_1s_mat)./size(all_sig_cells_1s_mat, 1);

    %plotting resp fracs across odors
    means_long = nanmean(saved_resp_vecs);
    means_1s = nanmean(saved_resp_vecs_1s);

    ses_long = nanstd(saved_resp_vecs)./sqrt(size(saved_resp_vecs, 1));
    ses_1s = nanstd(saved_resp_vecs_1s)./sqrt(size(saved_resp_vecs_1s, 1));

    
    %PLOTTING
    
    h_fig = figure(1)
    errorbar(means_long, ses_long, 'O', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 3)
    hold on
    errorbar(means_1s, ses_1s, 'O', 'Color', [0.7, 0.7, 0.7], 'LineWidth', 3)
    xlabel('Odor number')
    ylabel('Fraction of responsive cells')
    set(h_fig, 'units','normalized','position',[.1 .1 .25 .35])
    
    
    hold off

    all_resp_vecs_long = reshape(saved_resp_vecs, [], 1);
    all_resp_vecs_1s = reshape(saved_resp_vecs_1s, [], 1);
    
    axis([0, 6, 0, 1])

    [h, p] = ttest(all_resp_vecs_long, all_resp_vecs_1s);       %paired sample ttest
    disp(['Comparing responder fractions for 1s stim with longer stim fractions: p value ' num2str(p)])
    %------------------
    
    %checking if combination odor responses are randomly distributed
    
    %long odor responses
    %counting fractions of cells that respond to pairs of odors in real data
    pair_list = combinator(n_odors, 2, 'c');
    pair_resp_counts = zeros(size(pair_list, 1), 1); 
    for pair_n = 1:size(pair_list, 1)
        curr_pair = pair_list(pair_n, :);
        
        sum_vec = all_sig_cells_mat(:, curr_pair(1)) + all_sig_cells_mat(:, curr_pair(2));
        both_positives = find(sum_vec == 2);
        pair_resp_counts(pair_n, 1) = size(both_positives, 1);
    end
    pair_resp_counts = pair_resp_counts./size(all_sig_cells_mat, 1);
    
    [rand_dists, count_mat] = build_rand_dists(n_resp_vec, size(all_sig_cells_mat, 1), 1000);
    count_mat = count_mat./size(all_sig_cells_mat, 1);
    rand_means = mean(count_mat, 2);
    rand_sds = std(count_mat, [], 2);
    rand_sds = rand_sds.*2.33;
    
    
    fig_h = figure(6);
    shadedErrorBar([], rand_means, rand_sds, {'-', 'Color', [0.6, 0.6, 0.6]}, 1)
    hold on
    plot(pair_resp_counts, 'Ok', 'markerFaceColor', 'k')
    hold off
    ylabel('responder fractions')
    xlabel('odor-pairs')
    %title('Actual odor-pair responder fractions compared to those expected by random assignment.')
    set(fig_h, 'units','normalized','position',[.1 .1 .25 .35])
    %------------------
    
    
    
    %checking for fold increases in    
    h_fig = figure(2)    
    mean_ratios = mean(saved_resp_vec_ratio, 1);
    ses_ratios = std(saved_resp_vec_ratio)./sqrt(size(saved_resp_vec_ratio, 1));
    
    errorbar(mean_ratios, ses_ratios, 'O', 'Color', [0, 0, 0], 'LineWidth', 3)
    
    xlabel('Odor number')
    ylabel('Ratio of long stim responders to 1s stim responders')
    set(h_fig, 'units','normalized','position',[.1 .1 .25 .35])
        
    hold off

    axis([0, 6, 0, 1])
    
    keyboard
    
end

keyboard
%comparing responder fractions across cell-types
a = reshape(saved_resp_vecs_all, [], 1, n_direc_lists);
[p, tbl, stats] = anova1(a);
mult_comp_tbl = multcompare(stats);

a = reshape(saved_resp_vecs_1s_all, [], 1, n_direc_lists);
[p_1s, tbl, stats] = anova1(a);
mult_comp_tbl_1s = multcompare(stats);

