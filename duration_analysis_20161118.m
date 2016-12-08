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
    
    
    %parsing direc list path for name of direc list
    namei = findstr(list_direc, 'dataset_list');
    dir_list_name = (list_direc(namei:(end-4)));
        
    %loop to go through all experiment datasets listed in list file
    saved_traces_flies = [];
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
        
        long_dur_n = find(odor_dur_list == 60);
        
        for odor_n = 1:length(odor_list)
            odor_ni = odor_list(odor_n);
            
            %identifying significant responders to 60s stimuli, for current odor
            sig_cells = find(sig_cell_mat(:, odor_ni, long_dur_n) == 1);
            od_trs = find(stim_mat(:, 1) == odor_ni);
            
            stim_end_fr = stim_frame + ceil(60./frame_time);
            an_end_fr = min([n_frames, ceil(stim_end_fr + 200./frame_time)]);
            saved_traces_fly = zeros( n_frames, length(sig_cells), length(odor_dur_list) ) + nan;
            
            
            %building list of trs of current odor for each dur
            for dur_n = 1:length(odor_dur_list)
                dur_ni = odor_dur_list(dur_n);
                stim_end_fr = stim_frame + ceil(odor_dur_list(dur_n)./frame_time);
                an_end_fr = min([n_frames, ceil(stim_end_fr + 200./frame_time)]);
                dur_trs = find(stim_mat(:, 2) == dur_ni);
                curr_trs = intersect(od_trs, dur_trs);
                
                saved_traces_fly(:, :, dur_n) = nanmean(dff_data_mat(:, sig_cells, curr_trs, odor_ni), 3);
                
            end
            
            saved_traces_flies = concatenate_padded(saved_traces_flies, saved_traces_fly, 2, nan);      %pooling ave traces accross odors, flies. Keeping cells, od_dur distinct.
        
        end
    end
    
    %CLUSTERING
    %taking ave trace for cells pooled across odors and flies for 60s
    %stim and clustering
    response_matrix_r = squeeze(saved_traces_flies(stim_frame:an_end_fr, :, long_dur_n));
    n_cells = size(response_matrix_r, 2);

    %smoothing ave responses before clustering
    response_matrix = tsmovavg_m(response_matrix_r, 's', 10, 1);
    response_matrix(1:9, :) = response_matrix_r(1:9, :);            %replacing some baseline frames from un-smoothed matrix because tsmovavg returns nans there.
    
    %response_matrix(:, )
    un_clust = [];
    ccorr = corrcoef(response_matrix, 'rows', 'pairwise');
    
    save_path = 'C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Analysis\clustering_results\';
    if exist([save_path, dir_list_name, '.mat'], 'file') == 2
        clust_results = load([save_path, dir_list_name, '.mat']);
        clust_results = clust_results.clust_results;
        a = clust_results{1};
        b = clust_results{2};
        c = clust_results{3};
    
    elseif exist([save_path, dir_list_name, '.mat'], 'file') == 0
        %------------------------------------
        X = response_matrix';       %activity data matrix - cells x frames, with trials concatenated

        thr = 0.8;              %threshold of kmeans++ runs a pair of cells must co-occur in to be considered for clustering
        no_iters = 500;         %no of times meta_k_means_tank calls k_means_plusplus; set to 100 for preliminary corrcut scan
        clustat = [];
        un_clust = [];

        %loop to try out various values of corrcut
        ct_range = 0.65:0.05:0.9;
        for corrcut = 0.65:0.05:0.9
            %corrcut                             %corrcut is the threshold inter-cluster correlation coefficient above which they are merged
            [a b c] = meta_k_means_tank1(X, 'correlation', corrcut, thr, no_iters);


            clustcorr = [];
            clustnum = [];
            for i = 1 : size(a, 1)
                temp = a{i,1};
                clustnum = [clustnum length(temp)];
                tempcorr = 0;
                count = 1;
                for x = 1 : length(temp)-1
                    for y = x + 1 : length(temp)
                        tempcorr = tempcorr + ccorr(temp(x), temp(y));
                        count = count + 1;
                    end
                end
                clustcorr = [clustcorr, tempcorr/count];
            end
            clustat = [clustat; [mean(clustnum) (sum(clustcorr.*clustnum)/sum(clustnum))]];
        end

        %identifying best value of corrcut
        x = clustat(:, 1);
        x = (x - min(x))/(max(x)-min(x));
        y = clustat(:, 2);
        y = (y - min(y))/(max(y)-min(y));
        [null, corrcuti] = max(x.*y);

        corrcut = ct_range(corrcuti);

        clear null

        %running meta-k-means to identify final clusters with
        %well-chosen value of corrcut (threshold corrcoef to fuse meta-clusters)
        no_iters = 1000;            
        [a, b, c] = meta_k_means_tank1(X, 'correlation', corrcut, thr, no_iters);

        %saving results of clustering to file 
        clust_results = {a, b, c};
        save([save_path, dir_list_name, '.mat'], 'clust_results');
        
        beep

    else
    end
    
    
    
    %building list of cells and their group numbers
    no_clusts = size(a, 1);
    cell_gp_vec = zeros(n_cells, 1);
    for c_num = 1:no_clusts
        c_list = a{c_num, 1};
        cell_gp_vec(c_list, 1) = c_num;
    end

    clear no_clusts
    clear c_num
    clear c_list

    cell_gp_vec = [cell_gp_vec, (1:1:n_cells)'];
    cell_gp_vec_orig = cell_gp_vec;

    %sorting cells within groups as per their corrcoeffs with
    %group-averaged trace
    c_lengths = [];
    c_list_f = cell_gp_vec_orig(:, 1);                  %full list of all cells' cluster numbers
    


    response_matrix_o = [];
    %accounting for un-clustered cells
    c_list = find(c_list_f == 0);
    un_clust = [un_clust; length(c_list), n_cells];
    for c_no = 1:length(c_list)
        c_noi = c_list(c_no);

    end
    
    %sorting response matrix by corrcoef with cluster mean
    for clust_no = 1:max(cell_gp_vec_orig(:, 1));
        c_list = find(c_list_f == clust_no);            %cell numbers in old list that belong to current cluster
        %condition to skip clusters with < 5 cells in them
        if length(c_list) < 5
            for c_no = 1:length(c_list)
                c_noi = c_list(c_no);
            end

            continue
        else
        end                
        c_lengths = [c_lengths; length(c_list)]; 
        cr_list = c(c_list, clust_no);                  %list of corrcoeffs of cells in c_list, with their own cluster's averaged trace 
        curr_trace_mat = saved_traces_flies(:, c_list, long_dur_n)';       %matrix of activity data for each cell that belongs to this cluster

        curr_trace_mat = [cr_list, curr_trace_mat];     %concatenating corrcoefs to activity traces to sort both together
        curr_trace_mat = sortrows(curr_trace_mat, -1);      %sorting rows by first column ie corrcoefs
        response_matrix_o = [response_matrix_o; curr_trace_mat(:, 2:end)];
    end

    max_mat = repmat(nanmax(response_matrix_o, [], 2), 1, size(response_matrix_o, 2));
    response_matrix_o_norm = response_matrix_o./max_mat;
    response_matrix_o_norm_s = tsmovavg_m(response_matrix_o_norm, 's', 10, 2);
    response_matrix_o_norm_s(:, 1:9) = response_matrix_o_norm(:, 1:9);

    %calculating new corrcoeff mat with re-arranged cells
    curr_color = color_vec(2, :);
    figure(1)
    imagesc(response_matrix_o_norm_s, [0, 1]);
    stim_frs_saved = [stim_frame, stim_end_fr];
    add_stim_shading(1, stim_frs_saved, 0.20, curr_color)
    colormap(greymap)
    set_xlabels_time(1, frame_time, .4)
    ylabel('sig. cell-odor pairs')
    xlabel('time (s)')
    
    c_o = corrcoef(response_matrix_o', 'rows', 'pairwise');
    figure(2)
    imagesc(c_o)
    colormap('jet')
    xlabel('sig. cell-odor pairs')
    ylabel('sig. cell-odor pairs')
    
    figure(3)
    imagesc(ccorr)
    colormap('jet')
    xlabel('sig. cell-odor pairs')
    ylabel('sig. cell-odor pairs')
    %------------------------------------
    n_clusts = max(cell_gp_vec(:, 1));
    for clust_n = 0:n_clusts
        curr_clust_list = find(cell_gp_vec(:, 1) == (clust_n));
        curr_traces = squeeze(saved_traces_flies(:, curr_clust_list, long_dur_n));
        curr_traces = curr_traces./repmat(max(curr_traces), size(curr_traces, 1), 1);

        ave_resp_trace = nanmedian(curr_traces, 2); 

        figure(4)
        plot(curr_traces, 'LineWidth', 1, 'Color', [.75, .75, .75])
        hold on
        plot(ave_resp_trace, 'LineWidth', 2, 'Color', [0, .447, .741])
        hold off
        add_stim_shading(4, stim_frs_saved, 0.20, curr_color)

        xlabel('time')
        ylabel('Normalized dF/F')
        disp(['Cluster size ' num2str(length(curr_clust_list)./size(response_matrix, 2))])
        
        if clust_n == 0
            title(['un-clustered cells, n ' int2str(length(curr_clust_list)), '; ', int2str(length(curr_clust_list)./length(cell_gp_vec).*100), '%' ])
        else
            title(['n ' int2str(length(curr_clust_list)), '; ', int2str(length(curr_clust_list)./length(cell_gp_vec).*100), '%' ])
        end
        del = input('press enter');
        
    end


    ave_resp_trace = [];
    
    
    
%----------------------------------------------------------------
%CLASSIFYING CLUSTERS AS ONSET, SUSTAINED OR OFF RESPONSE CLUSTERS

    
    
    
    
    keyboard
end
    
        