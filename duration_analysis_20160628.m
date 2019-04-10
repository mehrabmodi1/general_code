clear all
close all

direc_lists_mat = [...%{'D:\Data\CSHL\dataset_list_sustained_MB418B_20160414.txt'};... %KC A'B'
                    %{'D:\Data\CSHL\dataset_list_sustained_MB185B_20160426.txt'};... %KC AB
                    %{'D:\Data\CSHL\dataset_list_stim_dur_20160316.txt'} ... KCs all
                    %{'D:\Data\CSHL\dataset_list_sustained_MB131B_20160515.txt'}... %KC G
                    {'D:\Data\CSHL\dataset_list_PN_GH146_20161002.txt'}... %PN axons
                    ]; 
                
       
            
n_direc_lists = size(direc_lists_mat, 1);
                
                
color_vec = load('C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Code\std_color_vec.txt');
a = colormap('bone');
greymap = flipud(a);

suppress_plots = 0;       %1 - doesn't plot stuff, 0 - plots stuff

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
    
    all_response_mat = [];
    all_sd_mat = [];
    all_ave_resp_mat = [];
    all_ave_sd_mat = [];
    all_saved_cell_vecs = [];
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
        
        clear raw_data_mat
        del = isnan(stim_mat(:, 1));
        del = find(del == 1);
        stim_mat(del, :) = [];
        dff_data_mat(:, :, del, :) = [];                    %getting rid of trials for which stim data is not available
        odor_list = unique(stim_mat(:, 1) );
        n_odors = length(odor_list);
        n_frames = size(dff_data_mat, 1);
        n_cells = size(dff_data_mat, 2);
        n_trials = size(stim_mat, 1);
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
        [resp_areas, sig_trace_mat, sig_cell_mat, sig_cell_block_mati, sig_cell_1s_mat] = ...
            cal_sig_responses_20160615(dataset, dff_data_mat, stim_mat, prot_switch_trials, list_direc, an_trial_window);

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
        
      
        %measuring response size for each cell with a small, moving window
        [response_mat, sd_mat, ave_resp_mat, ave_sd_mat] = find_pk_resp(dff_data_mat, stim_mat, frame_time, n_reps, stim_time);
        
        %concatenating cells across flies
        all_response_mat = pad_n_concatenate(all_response_mat, response_mat, 1, nan);
        all_sd_mat = pad_n_concatenate(all_sd_mat, sd_mat, 1, nan);
        all_ave_resp_mat = pad_n_concatenate(all_ave_resp_mat, ave_resp_mat, 1, nan);
        all_ave_sd_mat = pad_n_concatenate(all_ave_sd_mat, ave_sd_mat, 1, nan);
        
        
        %Analysing response time courses
        saved_cell_vecs = zeros(3, n_cells, n_odors) + nan;
        for cell_n = 1:n_cells
            for odor_n = 1:n_odors
                odor_dur_n = length(odor_dur_list);     %only looking at longest odor duration
                odor_dur = odor_dur_list(odor_dur_n);
                od_trs = find(stim_mat(:, 1) == odor_n);
                dur_trs = find(stim_mat(:, 2) == odor_dur);
                curr_trs = intersect(od_trs, dur_trs);
                start_frame = stim_frame;
                end_frame = round(odor_dur_list(odor_dur_n)./frame_time);
                curr_traces = squeeze(dff_data_mat(:, cell_n, curr_trs, odor_n));
                ave_tr = nanmean(curr_traces, 2);           %averaged dF/F trace for curr cell, curr od, curr od dur
                ave_tr = ave_tr(start_frame:end);
                ave_tr_f = tsmovavg_m(ave_tr, 's', 5, 1);
                
                %identifying location of peak and area after peak until
                %stim end
                [resp_size, pk_fr] = nanmax(ave_tr_f);
                plot(ave_tr_f)
                
                if pk_fr < end_frame
                    ave_area_post_n = nanmean(ave_tr_f(pk_fr:end_frame))./resp_size;    %normalised to size of pk resp
                    
                else
                    ave_area_post_n = 0;
                end
                
                cell_vec = [resp_size, pk_fr, ave_area_post_n];
                saved_cell_vecs(1:3, cell_n, odor_n) = cell_vec;
                
            end
            
        end
        all_saved_cell_vecs = pad_n_concatenate(all_saved_cell_vecs, saved_cell_vecs, 2, nan);
    end
    fclose(fid);
    
    del = find(all_ave_resp_mat < all_ave_sd_mat.*2.33);
    
    all_ave_resp_mat(del) = 0;                                      %forcing negative respnse values (less than noise floor) to 0          
    
    all_ave_resp_mat_norm = all_ave_resp_mat - repmat(all_ave_resp_mat(:, :, 1), [1, 1, size(all_ave_resp_mat, 3)]);        %subtracting 1s response from all responses
    
    
    %statistical testing
    %re-organising response size matrix so that responses for each odor are
    %treated as separate
    resh_ave_resp_mat = [];
    for odor_n = 1:n_odors
        resh_ave_resp_mat = pad_n_concatenate(resh_ave_resp_mat, squeeze(all_ave_resp_mat(:, odor_n, [1, 3])), 1, nan);
    end
    
       
    figure(1)
    fig_h = scattered_dot_plot(resh_ave_resp_mat, 1, 1, 5, 7, [0.5, 0.5, 0.5], 1, [0.8, 0.8, 0.8]);
    hold on
    
    means = nanmean(resh_ave_resp_mat);
    ses = nanstd(resh_ave_resp_mat)./sqrt(size(resh_ave_resp_mat, 1));
    
    errorbar([1, 7], means, ses, 'O', 'markerEdgeColor', [0.8, 0.4, 0.3], 'markerSize', 12, 'lineWidth', 4, 'Color', [0.8, 0.4, 0.3])
    hold off
    xlabel('stimulus duration (s)')
    ylabel('peak dF/F')
    ax = gca;
    ax.XTick = [1 7];
    ax.XTickLabel = {int2str(odor_dur_list(1)), int2str(odor_dur_list(3))};
    set(fig_h, 'units','normalized','position',[.1 .1 .25 .35])
    
    [p] = signrank(resh_ave_resp_mat(:, 1), resh_ave_resp_mat(:, 2));
    
    
    %plotting pop response distributions and calculating pop resp kurtosis
    %(sparseness measure from Willmore and Tolhurst, 2000)
    for odor_dur_n = 1:length(odor_dur_list)                
        hist_vecs = [];
        for odor_n = 1:n_odors        
            odor_dur_ni = odor_dur_n;
            
            all_response_mat_z = all_response_mat./all_sd_mat;      %resps z-scored to baseline sd-s
            curr_resps = all_response_mat_z(:, odor_n, odor_dur_ni, :);
            ave_resps = nanmean(squeeze(curr_resps), 2);        %averaging across repeats
            
            %calculating population kurtosis (sparseness measure)
            mean_resp = nanmean(ave_resps, 1);          %averaging across cells
            sd_resp = nanstd(ave_resps, 1);             %std across cells
            popz_ave_resps = (ave_resps - mean_resp)./sd_resp;
            pop_kurt = nanmean((popz_ave_resps.^4)) - 3;
            saved_kurts(odor_n, odor_dur_n, direc_list_n) = pop_kurt;
            
            
            curr_resps = reshape(curr_resps, [], 1);
            
            %plotting dist for curr odor
            figure(2)
            h = histogram(curr_resps, [0:0.5:10]);
            hist_vecs = [hist_vecs; h.Values];
        end
        
        figure(2)
        plot([0:0.5:9.5], hist_vecs)
        title(['odor duration ' int2str(odor_dur_list(odor_dur_ni))])
        xlabel('response size, z-scored to baseline')
        ylabel('number of responses')
       
    end
    
    
    keyboard
    
end


