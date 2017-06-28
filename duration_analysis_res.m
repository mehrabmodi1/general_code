clear all
close all

direc_lists_mat =  [{'D:\Data\CSHL\Resonant\dataset_list_mRFP_trip_trans_odor_durs.xls'};... %KC A'B'
                    {'D:\Data\CSHL\Resonant\dataset_list_mRFP_trip_trans_odor_pulses.xls'};... %KC AB
                   ]; 
                
            
n_direc_lists = size(direc_lists_mat, 1);
                
                
color_vec = load('C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Code\std_color_vec.txt');
a = colormap('bone');
greymap = flipud(a);

suppress_plots = 0;       %1 - doesn't plot stuff, 0 - plots stuff

saved_cell_data = cell(1, 45);

%figure property initialisation variables
plot_height = 200;
plot_width = 280;
axis_font_size = 8;
tick_length = [0.02, 0.05];
line_widths = 2;
marker_sizes = 1;
marker_sizes_f = 4;
box_l_width = 0.5;


saved_resp_vecs_all = [];
saved_resp_vecs_1s_all = [];

[del, odor_names] = xlsread('D:\Data\CSHL\odor_names_20161108.xls', 1);
n_tot_cells = 0;
n_tot_pairs = 0;
%loop to go through all directory lists
for direc_list_n = 1:n_direc_lists

    list_direc = direc_lists_mat{direc_list_n, 1};
    [del, curr_direc_list] = xlsread(list_direc, 1);
    n_dirs = size(curr_direc_list, 1);
    direc_counter = 0;
    
    %parsing direc list path for name of direc list
    namei = findstr(list_direc, 'dataset_list');
    dir_list_name = (list_direc(namei:(end-4)));
        
    %loop to go through all experiment datasets listed in list file
    for direc_counter = 1:n_dirs
        direc = curr_direc_list{direc_counter, 1};
        
        direc = [direc, '\'];
        
        %loading extracted raw fluorescence data matrices written by
        %raw_dff_extractor
        raw_data_mat = load([direc 'extracted_raw_data_mat.mat']);
        raw_data_mat = raw_data_mat.raw_data_mat;           %raw F traces extracted from ROIs
        
        
        %loading in and parsing params file to get stimulus paramater
        %details
        [stim_mat, stim_mat_simple, column_heads] = load_params_res(direc);
        odor_list = unique(stim_mat_simple(:, 2) );
        n_odors = length(odor_list);
        odor_dur_list = unique(stim_mat_simple(:, 3) );
        n_od_durs = length(odor_dur_list);
        
        %loading Suite2P results file
        Suite2P_results = load_suite2P_results(direc);
        frame_rate = Suite2P_results.ops.imageRate;
        frame_time = 1./frame_rate.*1000;     %in ms
        
        %calculating dF/F traces from raw data
        [dff_data_mat] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, direc);
        
        %identifying sig responses on a single trial basis, and then sig
        %responder cells for any given trial type
        [resp_areas, sig_trace_mat, sig_cell_mat, sig_trace_mat_old, sig_cell_mat_old] = ...
            cal_sig_responses_20161024(dataset, dff_data_mat, stim_mat, prot_switch_trials, list_direc, an_trial_window);
        
        
    end
end
        

%         
%         
%         %identifying sig responses on a single trial basis, and then sig
%         %responder cells in any individual block
%         [resp_areas, sig_trace_mat, sig_cell_mat, sig_trace_mat_old, sig_cell_mat_old] = ...
%             cal_sig_responses_20161024(dataset, dff_data_mat, stim_mat, prot_switch_trials, list_direc, an_trial_window);
%             
%         del = isnan(resp_areas(1, :));
%         bad_tr_list = find(del == 1);                   %list of trials thrown away due to movement
%         good_tr_list = 1:n_trials;
%         good_tr_list(bad_tr_list) = [];
%         n_prot_blocks = length(prot_switch_trials);     %number of protocol blocks in this dataset
% 
%         %making sure that at least half the trials in each block, for each odor were acquired
%         curr_window = n_trials./(length(odor_list).*n_prot_blocks);
%         bad_od_blk_counter = 0;
%         for odor_n = 1:length(odor_list)
%             odor_ni = odor_list(odor_n);
%             for odor_dur_n = 1:length(odor_dur_list)
%                 odor_dur_ni = odor_dur_list(odor_dur_n);
%                 for block_n = 1:n_prot_blocks
%                     od_blk_window_trs = odor_trial_list_builder_20160317(stim_mat, prot_switch_trials, odor_ni, odor_dur_ni, block_n, an_trial_window, 1);
% 
%                     n_bad_trs = length(intersect(od_blk_window_trs, bad_tr_list));            %number of bad trials in current block
%                     if n_bad_trs > curr_window/2
%                         disp('More than half the trials missing for one odor for one block.')
%                         bad_od_blk_counter = bad_od_blk_counter + 1;
%                     else
%                     end
% 
%                 end
%             end
%         end
%         if bad_od_blk_counter > 1
% 
%             %continue
%         else
%         end
%         
%         long_dur_n = find(odor_dur_list == 60);
%         
%         for odor_n = 1:length(odor_list)
%             odor_ni = odor_list(odor_n);
%             
%             %identifying significant responders to 60s stimuli, for current odor
%             sig_cells = find(sig_cell_mat(:, odor_ni, long_dur_n) == 1);
%             od_trs = find(stim_mat(:, 1) == odor_ni);
%             
%             stim_end_fr = stim_frame + ceil(60./frame_time);
%             an_end_fr = min([n_frames, ceil(stim_end_fr + 200./frame_time)]);
%             saved_traces_fly = zeros( n_frames, length(sig_cells), length(odor_dur_list) ) + nan;
%             saved_all_traces_fly = zeros( n_frames, n_cells, length(odor_dur_list) ) + nan;
%             
%             %building list of trs of current odor for each dur
%             for dur_n = 1:length(odor_dur_list)
%                 dur_ni = odor_dur_list(dur_n);
%                 stim_end_fr = stim_frame + ceil(odor_dur_list(dur_n)./frame_time);
%                 an_end_fr = min([n_frames, ceil(stim_end_fr + 200./frame_time)]);
%                 dur_trs = find(stim_mat(:, 2) == dur_ni);
%                 curr_trs = intersect(od_trs, dur_trs);
%                 
%                 saved_traces_fly(:, :, dur_n) = nanmean(dff_data_mat(:, sig_cells, curr_trs, odor_ni), 3);
%                 saved_all_traces_fly(:, :, dur_n) = nanmean(dff_data_mat(:, :, curr_trs, odor_ni), 3);
%             end
%             
%             
%             %plotting single cell traces for all sig cells one by one
%                 if suppress_plots == 0
%                     for sig_cell_n = 1:length(sig_cells)
%                         curr_cell_n = sig_cells(sig_cell_n);
%                         for dur_n = 1:length(odor_dur_list)
%                             dur_ni = odor_dur_list(dur_n);
%                             dur_trs = find(stim_mat(:, 2) == dur_ni);
%                             curr_trs = intersect(od_trs, dur_trs);
%                             stim_end_fr = stim_frame + ceil(dur_ni./frame_time);
%                             curr_traces_p = squeeze(dff_data_mat(:, curr_cell_n, curr_trs, odor_ni));
%                             ave_trace_p = nanmean(dff_data_mat(:, curr_cell_n, curr_trs, odor_ni), 3);
%                             norm_n = max(ave_trace_p);
%                             norm_n = 1;
%                             curr_traces_p = curr_traces_p./norm_n;
%                             ave_trace_p = ave_trace_p./norm_n;
%                             fig_h = figure(8 + dur_n);
%                             plot(curr_traces_p, 'Color', [.65, .65, .65])
%                             hold on
%                             plot(ave_trace_p, 'Color', 'k', 'LineWidth', 2)
%                             hold off
%                             ax = axis;
%                             if dur_n > 1
%                                 ax = max([ax; ax_old]);
%                             else
%                             end
%                             ax_old = ax;
%                             axis([0, ax(2), ax(3), 3.5]);
%                             ylabel('dF/F')
%                             set_xlabels_time((8 + dur_n), frame_time, .5)
%                             fig_wrapup((8+dur_n));
%                                                         
%                             
%                         end
%                         
%                         for dur_n = 1:length(odor_dur_list)
%                             figure((8 + dur_n))
%                             axis([0, ax(2), ax(3), ax(4)]);
%                             dur_ni = odor_dur_list(dur_n);
%                             stim_end_fr = stim_frame + ceil(dur_ni./frame_time);
%                             add_stim_bar((8+dur_n), [stim_frame, stim_end_fr], color_vec(odor_n, :));
% 
%                         end
%                         keyboard                        
%                         for dur_n = 1:length(odor_dur_list)
%                             close(figure(8 + dur_n))
%                         end
%                     end
%                 end
%             
%             saved_traces_flies = concatenate_padded(saved_traces_flies, saved_traces_fly, 2, nan);      %pooling ave traces accross odors, flies. Keeping cells, od_dur distinct.
%             saved_all_traces_flies = concatenate_padded(saved_all_traces_flies, saved_all_traces_fly, 2, nan);
%             saved_cell_ns = [saved_cell_ns; (sig_cells + prev_n_cells)];
%         end
%         n_tot_cells = [n_tot_cells + n_cells];
%         n_tot_pairs = [n_tot_pairs + (n_cells.*n_odors)];
%         prev_n_cells = prev_n_cells + n_cells;
%     end
%    
%     %CLUSTERING
%     %taking ave trace for cells pooled across odors and flies for 60s
%     %stim and clustering
%     response_matrix_r = squeeze(saved_traces_flies(stim_frame:an_end_fr, :, long_dur_n));
%     n_cells = size(response_matrix_r, 2);
% 
%     %smoothing ave responses before clustering
%     response_matrix = tsmovavg_m(response_matrix_r, 's', 10, 1);
%     response_matrix(1:9, :) = response_matrix_r(1:9, :);            %replacing some baseline frames from un-smoothed matrix because tsmovavg returns nans there.
%     
%     %response_matrix(:, )
%     un_clust = [];
%     ccorr = corrcoef(response_matrix, 'rows', 'pairwise');
%     
%     save_path = 'C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Analysis\clustering_results\';
%     if exist([save_path, dir_list_name, '.mat'], 'file') == 2
%         clust_results = load([save_path, dir_list_name, '.mat']);
%         clust_results = clust_results.clust_results;
%         a = clust_results{1};
%         b = clust_results{2};
%         c = clust_results{3};
%     
%     elseif exist([save_path, dir_list_name, '.mat'], 'file') == 0
%         %------------------------------------
%         X = response_matrix';       %activity data matrix - cells x frames, with trials concatenated
% 
%         thr = 0.8;              %threshold of kmeans++ runs a pair of cells must co-occur in to be considered for clustering
%         no_iters = 500;         %no of times meta_k_means_tank calls k_means_plusplus; set to 100 for preliminary corrcut scan
%         clustat = [];
%         un_clust = [];
% 
%         %loop to try out various values of corrcut
%         ct_range = 0.65:0.05:0.9;
%         for corrcut = 0.65:0.05:0.9
%             %corrcut                             %corrcut is the threshold inter-cluster correlation coefficient above which they are merged
%             [a b c] = meta_k_means_tank1(X, 'correlation', corrcut, thr, no_iters);
% 
% 
%             clustcorr = [];
%             clustnum = [];
%             for i = 1 : size(a, 1)
%                 temp = a{i,1};
%                 clustnum = [clustnum length(temp)];
%                 tempcorr = 0;
%                 count = 1;
%                 for x = 1 : length(temp)-1
%                     for y = x + 1 : length(temp)
%                         tempcorr = tempcorr + ccorr(temp(x), temp(y));
%                         count = count + 1;
%                     end
%                 end
%                 clustcorr = [clustcorr, tempcorr/count];
%             end
%             clustat = [clustat; [mean(clustnum) (sum(clustcorr.*clustnum)/sum(clustnum))]];
%         end
% 
%         %identifying best value of corrcut
%         x = clustat(:, 1);
%         x = (x - min(x))/(max(x)-min(x));
%         y = clustat(:, 2);
%         y = (y - min(y))/(max(y)-min(y));
%         [null, corrcuti] = max(x.*y);
% 
%         corrcut = ct_range(corrcuti);
% 
%         clear null
% 
%         %running meta-k-means to identify final clusters with
%         %well-chosen value of corrcut (threshold corrcoef to fuse meta-clusters)
%         no_iters = 1000;            
%         [a, b, c] = meta_k_means_tank1(X, 'correlation', corrcut, thr, no_iters);
% 
%         %saving results of clustering to file 
%         clust_results = {a, b, c};
%         save([save_path, dir_list_name, '.mat'], 'clust_results');
%         
%         beep
% 
%     else
%     end
%     
%     
%     
%     %building list of cells and their group numbers
%     no_clusts = size(a, 1);
%     cell_gp_vec = zeros(n_cells, 1);
%     for c_num = 1:no_clusts
%         c_list = a{c_num, 1};
%         cell_gp_vec(c_list, 1) = c_num;
%     end
% 
%     clear no_clusts
%     clear c_num
%     clear c_list
% 
%     cell_gp_vec = [cell_gp_vec, (1:1:n_cells)'];
%     cell_gp_vec_orig = cell_gp_vec;
% 
%     %sorting cells within groups as per their corrcoeffs with
%     %group-averaged trace
%     c_lengths = [];
%     c_list_f = cell_gp_vec_orig(:, 1);                  %full list of all cells' cluster numbers
%     
% 
% 
%     response_matrix_o = [];
%     response_matrix_o_20 = [];
%     %accounting for un-clustered cells
%     c_list = find(c_list_f == 0);
%     un_clust = [un_clust; length(c_list), n_cells];
%     for c_no = 1:length(c_list)
%         c_noi = c_list(c_no);
% 
%     end
%     
%     
%     %computing time of center of mass of ave response and classifying
%     %cluster as onset or sus with a CoM cut-off time of 60/3 = 20s
%     frame_vec = 1:1:(stim_end_fr - stim_frame);
%     saved_CoMs = [];
%     for clust_n = 1:max(cell_gp_vec)
%         curr_clust_list = find(cell_gp_vec(:, 1) == (clust_n));
%         curr_traces = squeeze(saved_traces_flies(:, curr_clust_list, long_dur_n));
%         curr_traces = curr_traces./repmat(max(curr_traces), size(curr_traces, 1), 1);
%         curr_ave = nanmedian(curr_traces, 2);
%         ave_resp_traces(:, (clust_n + 1)) = curr_ave;
%         
%         %computing CoM of current ave trace
%         curr_avei = curr_ave((stim_frame + 1):stim_end_fr);
%         del = find(curr_avei < .4);
%         curr_avei(del) = 0;
%         CoM = nansum(curr_avei.*frame_vec')./nansum(curr_avei);         %just like calculating the mean of a frequency distribution
%         if CoM > (20./frame_time)
%             sus_c = 1;
%         else
%             sus_c = 0;
%         end
%         
%         saved_CoMs = [saved_CoMs; CoM, clust_n, sus_c];
%     end
%     
%     
%     %sorting response matrix by corrcoef with cluster mean
%     for clust_no = 1:max(cell_gp_vec_orig(:, 1));
%         c_list = find(c_list_f == clust_no);            %cell numbers in old list that belong to current cluster
%                 
%         %condition to skip clusters with < 5 cells in them
%         if length(c_list) < 5
%             for c_no = 1:length(c_list)
%                 c_noi = c_list(c_no);
%             end
%             
% 
%             continue
%         else
%         end                
%         c_lengths = [c_lengths; length(c_list)]; 
%         cr_list = c(c_list, clust_no);                  %list of corrcoeffs of cells in c_list, with their own cluster's averaged trace 
%         curr_trace_mat = saved_traces_flies(:, c_list, long_dur_n)';       %matrix of activity data for each cell that belongs to this cluster
%         curr_trace_mat_20 = saved_traces_flies(:, c_list, (long_dur_n-1) )';             %matrix of activity data for 20s odor trials 
%         
%         curr_trace_mat = [cr_list, curr_trace_mat];     %concatenating corrcoefs to activity traces to sort both together
%         curr_trace_mat_20 = [cr_list, curr_trace_mat_20];
%         
%         curr_trace_mat = sortrows(curr_trace_mat, -1);      %sorting rows by first column ie corrcoefs
%         curr_trace_mat_20 = sortrows(curr_trace_mat_20, -1);
%         
%         response_matrix_o = [response_matrix_o; curr_trace_mat(:, 2:end)];
%         response_matrix_o_20 = [response_matrix_o_20; curr_trace_mat_20(:, 2:end)];
%     end
% 
%     max_mat = repmat(nanmax(response_matrix_o, [], 2), 1, size(response_matrix_o, 2));
%     response_matrix_o_norm = response_matrix_o./max_mat;
%     
%     max_mat = repmat(nanmax(response_matrix_o_20, [], 2), 1, size(response_matrix_o_20, 2));
%     response_matrix_o_norm_20 = response_matrix_o_20./max_mat;
%     
%     if isempty(response_matrix_o_norm) == 1
%         continue
%     else
%     end
%     response_matrix_o_norm_s = tsmovavg_m(response_matrix_o_norm, 's', 10, 2);
%     response_matrix_o_norm_s(:, 1:9) = response_matrix_o_norm(:, 1:9);
% 
%     response_matrix_o_norm_20_s = tsmovavg_m(response_matrix_o_norm_20, 's', 10, 2);
%     response_matrix_o_norm_20_s(:, 1:9) = response_matrix_o_norm_20(:, 1:9);
%     
%     %calculating new corrcoeff mat with re-arranged cells
%     curr_color = color_vec(2, :);
%     fig_h = figure(1);
%     imagesc(response_matrix_o_norm_s, [0, 1]);
%     stim_frs_saved = [stim_frame, stim_end_fr];
%     colormap(greymap)
%     set_xlabels_time(1, frame_time, .5)
%     ylabel('sig. cell-odor pairs')
%     xlabel('time (s)')
%     set(fig_h, 'Position', [100, 100, 100 + plot_width, 100 + plot_height]);
%     fig_wrapup(1);
%     add_stim_bar(1, stim_frs_saved, curr_color)
%     
%     
%     fig_h = figure(11);
%     n_real_frames = n_frames - (sum(isnan(response_matrix_o_norm_20_s(1, :))));
%     imagesc(response_matrix_o_norm_20_s(:, 9:n_real_frames), [0, 1]);
%     stim_frs_short = [(stim_frame - 9), round((stim_frame - 9) + odor_dur_list((long_dur_n-1))./frame_time)];
%     colormap(greymap)
%     set_xlabels_time(11, frame_time, .51)
%     ylabel('sig. cell-odor pairs')
%     xlabel('time (s)')
%     fig_wrapup(11);
%     set(fig_h, 'Position', [100, 100, 100 + (plot_width.*(n_real_frames./n_frames)), 100 + plot_height]);
%     add_stim_bar(11, stim_frs_short, color_vec(2, :));
%     
%     
%     c_o = corrcoef(response_matrix_o', 'rows', 'pairwise');
%     fig_h = figure(2);
%     imagesc(c_o)
%     colormap('jet')
%     xlabel('sig. cell-odor pairs')
%     ylabel('sig. cell-odor pairs')
%     set(fig_h, 'Position', [100, 100, 100 + plot_width, 100 + plot_height]);
%     fig_wrapup(2);
%     
%     fig_h = figure(3);
%     imagesc(ccorr)
%     colormap('jet')
%     xlabel('sig. cell-odor pairs')
%     ylabel('sig. cell-odor pairs')
%     set(fig_h, 'Position', [100, 100, 100 + plot_width, 100 + plot_height]);
%     fig_wrapup(3);
%     
%     %------------------------------------
%     n_clusts = max(cell_gp_vec(:, 1));
%     for clust_n = 0:n_clusts
%         curr_clust_list = find(cell_gp_vec(:, 1) == (clust_n));
%         curr_traces = squeeze(saved_traces_flies(:, curr_clust_list, long_dur_n));
%         curr_traces = curr_traces./repmat(max(curr_traces), size(curr_traces, 1), 1);
% 
%         ave_resp_traces(:, (clust_n + 1)) = nanmedian(curr_traces, 2); 
% 
%         fig_h = figure(4);
%         plot(curr_traces, 'LineWidth', 1, 'Color', [.75, .75, .75])
%         hold on
%         plot(ave_resp_traces(:, (clust_n + 1)), 'LineWidth', 2, 'Color', [0, .447, .741])
%         hold off
%         add_stim_shading(4, stim_frs_saved, 0.20, curr_color)
%         set_xlabels_time(4, frame_time, .5)
%         fig_wrapup(4);
%         
%         xlabel('time')
%         ylabel('normalized dF/F')
%         %disp(['Cluster size ' num2str(length(curr_clust_list)./size(response_matrix, 2))])
%         set(fig_h, 'Position', [100, 100, 100 + plot_width, 100 + plot_height]);
%         
%         
%         if clust_n == 0
%             title(['un-clustered cells, n ' int2str(length(curr_clust_list)), '; ', int2str(length(curr_clust_list)./length(cell_gp_vec).*100), '%' ])
%         else
%             title(['n ' int2str(length(curr_clust_list)), '; ', int2str(length(curr_clust_list)./length(cell_gp_vec).*100), '%, ' int2str((length(curr_clust_list)./(n_tot_cells.*5)).*100) '% of all cell-odor pairs'])
%         end
%         
%         if suppress_plots == 0
%             del = input('press enter');
%             %keyboard
%         else
%         end
%     end
%     
%     ave_resp_traces(:, 1) = [];         %getting rid of un-clustered cells' trace
%     
%     %----------------------------------------------------------------
%     %Not using this criterion anymore
%     %CLASSIFYING CLUSTERS AS ONSET, SUSTAINED, OFF OR RAMP RESPONSE CLUSTERS
%     %calculating time to peak, area ratio for each smoothed ave resp trace
%     traces = saved_traces_flies(:, :, long_dur_n);
%     maxmat = repmat(nanmax(traces, [], 1), size(traces, 1), 1);
%     traces = traces./maxmat;                            %normalising
%     traces_s = tsmovavg_m(traces, 's', 10, 1);
%     traces_s(1:9, :) = traces(1:9, :);
%     traces_s = traces_s';
%     
%     [del pk_times] = nanmax(traces_s, [], 2);
%     pk_times = pk_times - stim_frame;
%     pk_times = pk_times.*frame_time;        %times to peak in s
%     
%     bin_edges = [0:6:60, inf];
%     bin_centers = bin_edges(1:11) + 3;
%     pk_dist_all = histcounts(pk_times, bin_edges);
%     pk_dist_all = pk_dist_all./sum(pk_dist_all).*100;
%     
%     
%     mid_frame = floor((stim_end_fr - stim_frame)./2) + stim_frame;
%     early_areas = nanmax(traces_s(:, stim_frame:mid_frame), [], 2);
%     early_areas = MinMaxCheck(repmat(0, 1, length(early_areas)), repmat(Inf, 1, length(early_areas)), early_areas);
%     late_areas = nanmax(traces_s(:, (mid_frame + 1):stim_end_fr), [], 2);
%     late_areas = MinMaxCheck(repmat(0, 1, length(late_areas)), repmat(Inf, 1, length(late_areas)), late_areas);
%     area_ratios = late_areas./early_areas;
%     
%     saved_dists = zeros((length(bin_edges) - 1), n_clusts);
%     %calculating overlap of pk time distributions and fusing clusters with
%     %more than 33.33% overlap
%     for clust_n = 1:n_clusts
%         curr_cells = find(c_list_f == clust_n);
%         curr_dist = histcounts(pk_times(curr_cells), bin_edges);
%         curr_dist = curr_dist./sum(curr_dist).*100;
%         saved_dists(:, clust_n) = curr_dist; 
%     end
%     
%     ovlap_pairs = [];
%     for clust_n = 1:(n_clusts - 1)
%         for clust_ni = (clust_n + 1):n_clusts
%             [overlap_areas] = dist_overlap_calc(bin_centers, saved_dists(:, clust_n), saved_dists(:, clust_ni));
%             overlap_area = sum(overlap_areas);
%             
%             if overlap_area >  40
%                 %checking if either member of this clust pair has already been pooled with another cluster
%                 ovlap_pairs = [ovlap_pairs; clust_n, clust_ni];
%             end
%             
%         end
%     end
%     
%    
%     %fusing all clusters that have more than 33.33% dist area overlap
%     saved_pools = [];
%     for clust_n = 1:n_clusts
%         saved_pools{(size(saved_pools, 2) + 1)} = clust_n;
%     end
%     for pair_n = 1:size(ovlap_pairs, 1)
%         saved_pools{(size(saved_pools, 2) + 1)} = ovlap_pairs(pair_n, :);
%     end
%     saved_pools_orig = saved_pools;
%     %checking if any of the newly created pools have any overlap and fusing them
%     still_fusing = 1;
%     counter = 1;
%     while still_fusing == 1
%         n_pools = size(saved_pools, 2);
%         new_pools = 0;
%         if n_pools > 1
%             saved_pools_new = saved_pools;
%             break_loops = 0;
%             for pool_n = 1:(n_pools - 1)
%                 curr_pool = saved_pools{pool_n};
%                 for pool_ni = (pool_n + 1):n_pools
%                     curr_pooli = saved_pools{pool_ni};
%                     if isempty(intersect(curr_pool, curr_pooli)) == 0
%                         saved_pools_new{(size(saved_pools_new, 2) + 1)} = unique(union(curr_pool, curr_pooli));
%                         ex_list = 1:1:size(saved_pools_new, 2);
%                         rem_list = [pool_n, pool_ni];
%                         ex_list(rem_list) = [];
%                         saved_pools_new = saved_pools_new(1, ex_list);
%                         new_pools = new_pools + 1;
%                         break_loops = 1;
%                         break                        
%                     else
%                     end
%                     
%                 end
%                 if break_loops == 1
%                     break
%                 else
%                 end
%             end
%             if new_pools > 0
%                 saved_pools = saved_pools_new;
%                 
%             elseif new_pools == 0
%                 still_fusing = 0;
%             else
%             end
%         else
%             still_fusing = 0;
%         end
%     end
%     
%     %changing cluster number to pool number for each cell acc to which pool
%     %it belongs
%     c_list_f_old = c_list_f;
%     c_list_f_new = c_list_f;
%     for pool_n = 1:size(saved_pools, 2)
%         curr_pool = saved_pools{pool_n};
%         n_clusts = size(curr_pool, 2);      %number of clusters in current pool
%         
%         for clust_n = 1:n_clusts
%             clust_ni = curr_pool(clust_n);
%             curr_cells = find(c_list_f == clust_ni);
%             c_list_f_new(curr_cells) = pool_n;
%         end
%         
%     end
%     c_list_f = c_list_f_new;
%     n_clusts = max(c_list_f);
%     %------------------- Done Fusing
%     
%     
%     
%     
%     %plotting
%     fig_h5 = figure(5);
%     plot(area_ratios, pk_times, '.k', 'MarkerSize', 30)
%     hold on
%     
%     fig_h6 = figure(6);
%     plot(bin_centers, pk_dist_all, 'k', 'LineWidth', 2)
%     hold on
%     
%     for clust_n = 1:n_clusts
%         curr_cells = find(c_list_f == clust_n);
%         curr_color = color_vec(clust_n, :);
% 
%         figure(5)
%         plot(area_ratios(curr_cells), pk_times(curr_cells), '.', 'Color', curr_color, 'MarkerSize', 30)
% 
%         figure(6)
%         plot(bin_centers, saved_dists(:, clust_n), 'Color', curr_color, 'LineWidth', 2)
%     end
%     figure(5)
%     fig_wrapup(5);
%     hold off
%     figure(6)
%     fig_wrapup(6);
%     hold off
%     
%     set(fig_h5, 'Position', [100, 100, 100 + plot_width, 100 + plot_height]);
%     set(fig_h6, 'Position', [100, 100, 100 + plot_width, 100 + plot_height]);
%     
%         
%     fig_h7 = figure(7);
%     %calculating grand average 60s response trace across all cells, across
%     %all odors
%     mean_trace = nanmean(squeeze(saved_all_traces_flies(:, :, long_dur_n)), 2);
%     ses = nanstd(squeeze(saved_all_traces_flies(:, :, long_dur_n)), [], 2)./sqrt(size(saved_all_traces_flies, 2));
%     time_vec = 0:frame_time:((length(mean_trace)-1).*frame_time);
%     shadedErrorBar(time_vec, mean_trace, ses, {'Color', [.9, .4, .5]});
%     xlabel('time (s)')
%     ylabel('dF/F')
%     add_stim_shading(7, [round(stim_frame.*frame_time), round(stim_end_fr.*frame_time)], 0.2, color_vec())
%     set(fig_h7, 'Position', [100, 100, 100 + plot_width, 100 + plot_height]);
%     fig_wrapup(7);
%     
%     %Counting numbers of onset cells, sus cells and cells that are both
%     n_cells_total = prev_n_cells;                           %total number of cells imaged i.e. list of all cell numbers (not just sig responders)
%     cell_type_vec = zeros(n_cells_total, 2) + nan;
%     for sus_type = 0:1
%         curr_clusts = find(saved_CoMs(:, 3) == sus_type);   %clust numbers identified as onset (0) or sus (1)
% 
%         %making list of all trace numbers that belong to any of the clusters in
%         %curr_clusts
%         all_trace_ns_curr = [];
%         for c_clustn = 1:length(curr_clusts)
%             c_clustni = curr_clusts(c_clustn);
%             curr_cs = find(c_list_f == c_clustni);
%             all_trace_ns_curr = [all_trace_ns_curr; curr_cs];
%         end
%         all_cells_curr = saved_cell_ns(all_trace_ns_curr);
%         all_cells_curr = unique(all_cells_curr);
% 
%         cell_type_vec(all_cells_curr, (sus_type + 1)) = 1;
%         
%     end
%     n_cells_onset = nansum(cell_type_vec(:, 1))./n_cells_total;
%     n_cells_sus = nansum(cell_type_vec(:, 2))./n_cells_total;
%     both_vec = nansum(cell_type_vec, 2);
%     n_cells_both = length(find(both_vec == 2))./n_cells_total;
%     n_cells_both_exp = n_cells_onset.*n_cells_sus;
%     fig_h8 = figure(8);
%     venn([n_cells_onset, n_cells_sus, 1], [n_cells_both, n_cells_onset, n_cells_sus, n_cells_both], 'FaceColor',{color_vec(1, :),color_vec(2, :), [1, 1, 1]},'FaceAlpha',{.5,.5,.5},'EdgeColor',{[.75, .75, .75], [.75, .75, .75], [.75, .75, .75]});
%     set(fig_h8, 'Position', [100, 100, 100 + plot_width, 100 + plot_height]);
%     fig_wrapup(8);
%     
%     %hold on
%     
%     
%     keyboard
%     ave_resp_traces = [];

    


        