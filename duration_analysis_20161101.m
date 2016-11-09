clear all
close all

direc_lists_mat =  [{'D:\Data\CSHL\dataset_list_sustained_MB418B_20160414.txt'};... %KC A'B'
                    {'D:\Data\CSHL\dataset_list_sustained_MB185B_20160426.txt'};... %KC AB
                    %{'D:\Data\CSHL\dataset_list_stim_dur_20160316.txt'}; ... KCs all
                    {'D:\Data\CSHL\dataset_list_sustained_MB131B_20160515.txt'};... %KC G
                    {'D:\Data\CSHL\dataset_list_stim_dur_20160316.txt'}... %OK107 KCs
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

[del, odor_names] = xlsread('D:\Data\CSHL\odor_names_20161108.xls', 1);

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
        
        direc = fgetl(fid);

        if ischar(direc) ~= 1
            break
        else
        end
        
        direc_counter = direc_counter + 1;

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
        
%         dump_direc = dump_direcs{direc_list_n};
%         save([dump_direc 'dFF_data_fly' int2str(direc_counter) '.mat'], 'dff_data_mat');
%         save([dump_direc 'stim_info_fly' int2str(direc_counter) '.mat'], 'stim_mat');
        
        clear raw_data_mat
        odor_list = unique(stim_mat(:, 1) );
        del = find(odor_list == 0);
        odor_list(del) = [];
        del = find(isnan(odor_list) == 1);
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
        del = find(odor_dur_list == 0);
        odor_dur_list(del) = [];
        del = find(isnan(odor_dur_list) == 1);
        odor_dur_list(del) = [];
        n_reps = dataset(1).stim.reps;
        stim_end_frames = [round(odor_dur_list./frame_time) + stim_frame];
        
        
        %identifying sig responses on a single trial basis, and then sig
        %responder cells in any individual block
        [resp_areas, sig_trace_mat, sig_cell_mat, sig_trace_mat_old, sig_cell_mat_old] = ...
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
        
%         %Plotting sig cell mats

%         figure(1)
%         subplot(3, 2, 1)
%         imagesc(sig_cell_mat(:, :, 1))
%         title('1s')
%         
%         subplot(3, 2, 3)
%         imagesc(sig_cell_mat(:, :, 2))
%         title('20s')
%         
%         subplot(3, 2, 5)
%         imagesc(sig_cell_mat(:, :, 3))
%         title('60s')
%         
%         %plotting old criterion stuff
%         subplot(3, 2, 2)
%         imagesc(sig_cell_mat_old(:, :, 1))
%         title('1s old')
%         
%         subplot(3, 2, 4)
%         imagesc(sig_cell_mat_old(:, :, 2))
%         title('20s old')
%         
%         subplot(3, 2, 6)
%         imagesc(sig_cell_mat_old(:, :, 3))
%         title('60s old')
        

        %plotting summary statistics and doing statistical tests on sig cell
        %frequencies
        for dur_n = 1:2
            if dur_n == 1
                odor_dur_n = find(odor_dur_list == 1);
            elseif dur_n == 2
                odor_dur_n = find(odor_dur_list == 60);
            else
            end
            
            try
                cell_freq_lists(:, odor_dur_n, direc_counter) = nanmean(sig_cell_mat(:, odor_list, odor_dur_n), 1);
                
            catch
                keyboard
            end
            
            
            
            
        end
        %keyboard
     
        
        %computing mean response area for each cell across repeats for 1s and
        %60s trials
        ave_areas = zeros(n_cells, length(odor_list), length(odor_dur_list)) + nan;
        se_areas = ave_areas;
        for odor_n = 1:length(odor_list)
            odor_ni = odor_list(odor_n);
            curr_odor_trs = find(stim_mat(:, 1) == odor_ni);
            for dur_n = 1:length(odor_dur_list)
                dur_ni = odor_dur_list(dur_n);
                curr_dur_trs = find(stim_mat(:, 2) == dur_ni);   
                curr_trs = intersect(curr_odor_trs, curr_dur_trs);
                ave_areas(:, odor_n, dur_n) = nanmean(resp_areas(:, curr_trs), 2); 
                se_areas(:, odor_n, dur_n) = nanstd(resp_areas(:, curr_trs), [], 2)./sqrt(length(curr_trs));

            end
        end
        dur_n = find(odor_dur_list == 1);
        ave_areas_1s = reshape(ave_areas(:, :, dur_n), 1, []);
        ses_areas_1s = reshape(se_areas(:, :, dur_n), 1, []);
        sig_cell_list_1s = find(sum(sig_cell_mat(:, :, dur_n), 2) > 0);
        dur_n = find(odor_dur_list == 60);
        ave_areas_60s = reshape(ave_areas(:, :, dur_n), 1, []);
        ses_areas_60s = reshape(se_areas(:, :, dur_n), 1, []);
        sig_cell_list_60s = find(sum(sig_cell_mat(:, :, dur_n), 2) > 0); 
        sig_cell_list = union(sig_cell_list_1s, sig_cell_list_60s);
        ave_areas = ave_areas(sig_cell_list, :);
        
        saved_ave_areas = concatenate_mat(saved_ave_areas, ave_areas, 1);


        
    end
    
    %COMPUTATIONS
    %1s mean, se cell-fracs
    dur_n = find(odor_dur_list == 1);
    means_1s = nanmean(cell_freq_lists(:, dur_n, :), 3);
    ses_1s = nanstd(cell_freq_lists(:, dur_n, :), [], 3)./sqrt(size(cell_freq_lists, 3));
    
    %60s mean, se cell-fracs
    dur_n = find(odor_dur_list == 60);
    means_60s = nanmean(cell_freq_lists(:, dur_n, :), 3);
    ses_60s = nanstd(cell_freq_lists(:, dur_n, :), [], 3)./sqrt(size(cell_freq_lists, 3));
    
    %computing change factor
    dur_n = find(odor_dur_list == 1);
    cell_freqs_1s = squeeze(cell_freq_lists(:, dur_n, :));
    dur_n = find(odor_dur_list == 60);
    cell_freqs_60s = squeeze(cell_freq_lists(:, dur_n, :));
    c_facs = (cell_freqs_60s - cell_freqs_1s)./(cell_freqs_60s + cell_freqs_1s);
    mean_facs = nanmean(c_facs, 2);
    ses_facs = nanstd(c_facs, [], 2)./sqrt(size(c_facs, 2));
    
   
    
    
    
    clear cell_freq_lists
    %PLOTTING
    %plotting responder fractions for each odor, for 1s and 60s stimuli
    figure(1)
    h_fig = figure(1);
    errorbar(means_60s, ses_60s, 'O', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 3)
    hold on
    errorbar(means_1s, ses_1s, 'O', 'Color', [0.7, 0.7, 0.7], 'LineWidth', 3)
    xticks = {odor_names{odor_list}};
    ax = gca;
    ax.XTick = 1:length(odor_list);
    ax.XTickLabel = xticks;
    ax.XTickLabelRotation=45;
    ylabel('Fraction of responsive cells')
    set(h_fig, 'units','normalized','position',[.1 .1 .25 .35])
    hold off
   
    %plotting change index (x-y)/(x+y)
    figure(2)
    h_fig = figure(2);
    errorbar(mean_facs, ses_facs, 'O', 'Color', 'k', 'LineWidth', 3)
    xticks = {odor_names{odor_list}};
    ax = gca;
    ax.XTick = 1:length(odor_list);
    ax.XTickLabel = xticks;
    ax.XTickLabelRotation=45;
    ylabel('Fraction change index')
    set(h_fig, 'units','normalized','position',[.1 .1 .25 .35])
    axis([0, (length(odor_list) + 1), -1, 1])
    hold off
    
    %plotting resp areas for 1s vs 60s stimuli
    figure(3)
    h_fig = figure(3);
    errorbar(ave_areas_1s(sig_cell_list), ave_areas_60s(sig_cell_list), ses_areas_60s(sig_cell_list), 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 4)
    hold on
    herrorbar(ave_areas_1s(sig_cell_list), ave_areas_60s(sig_cell_list), ses_areas_1s(sig_cell_list), '.k')
    axis_old = axis;
    min_ax = min([axis_old(1), axis_old(3)]);
    max_ax = max([axis_old(2), axis_old(4)]);
    diag_vec = [min_ax, max_ax];
    plot(diag_vec, diag_vec, '--', 'Color', '[0.75, 0.75, 0.75]', 'LineWidth', 1)
    hold off
    axis([min_ax, max_ax, min_ax, max_ax]);
    xlabel('1s odor responses (dF/F)');
    ylabel('60s odor responses (dF/F)');
    set(h_fig, 'units','normalized','position',[.1 .1 .25 .35])
    
    
    keyboard
    
end
    
        