clear all
close all

list_direc = 'D:\Data\CSHL\dataset_list_stim_dur_20160316.txt'; %expt datasets
color_vec = load('C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Code\std_color_vec.txt');
a = colormap('bone');
greymap = flipud(a);
colormap(greymap)


suppress_plots = 01;       %1 - doesn't plot stuff, 0 - plots stuff
fid = fopen(list_direc);
direc_counter = 0;

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

cell_counter = 0;
bigger_matrix = [];
%loop to go through all experiment datasets listed in list file
while 1
    direc_counter = direc_counter + 1;
    direc = fgetl(fid);
    
    if ischar(direc) ~= 1
        break
    else
    end
    
    direc = [direc, '\'];
    
    %replacing C: with D:
    a = strfind(direc, 'C:');
    direc(a) = 'D';
    
    
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
    odor_list = unique(stim_mat(:, 1) );
    n_odors = length(odor_list);
    n_frames = size(dff_data_mat, 1);
    n_cells = size(dff_data_mat, 2);
    n_trials = size(dff_data_mat, 3);
    stim_time = dataset(1).stim.stimLatency + 0.54;     %odor stim delay measured with PID is 540 ms
    frame_time = dataset(1).info.framePeriod;
    stim_frame = floor(stim_time./frame_time);
    an_trial_window = nan;
    odor_dur_list = unique(stim_mat(:, 2) );
    
    
    %identifying sig responses on a single trial basis, and then sig
    %responder cells in any individual block
    [resp_areas, sig_trace_mat, sig_cell_mat, sig_cell_block_mati] = ...
        cal_sig_responses_20160404(dataset, dff_data_mat, stim_mat, prot_switch_trials, list_direc, an_trial_window);
        
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
    
    
    %plotting ave response traces of all cells
    for odor_n = 1:length(odor_list)
       odor_ni = odor_list(odor_n);
       for odor_dur_n = 1:length(odor_dur_list)
           odor_dur_ni = odor_dur_list(odor_dur_n);
           curr_color = color_vec(odor_ni, :);
           stim_end_frm = ceil(stim_frame + odor_dur_ni./frame_time);
           %----------------- loops to walk through each odor, delivered for each duration
           
           tr_list = odor_trial_list_builder_20160317(stim_mat, prot_switch_trials, odor_ni, odor_dur_ni, block_n, an_trial_window, 1);               %list of trials of this odor, given for this duration
           ave_dff_mat = squeeze(dff_data_mat(:, :, tr_list, odor_ni) );
           ave_dff_mat = nanmean(ave_dff_mat, 3);          
           
           if suppress_plots == 0
               figure(1)
               imagesc(squeeze(ave_dff_mat)', [0, 4])
               xlabel('frame number')
               ylabel('cell number')
               set(gcf, 'Color', 'w')
               title(['fly number ' int2str(direc_counter) ', odor number ' int2str(odor_ni) ', odor duration' int2str(odor_dur_list(odor_dur_n))]);
               add_stim_shading(1, [stim_frame, stim_end_frm], 0.20, curr_color)
               set_xlabels_time(1, frame_time, 0.5)
               
               del = input('press enter');
           else
           end
       end
    end
   
    clear big_matrix
    %loop to analyse each cell across odors and durations
    for cell_n = 1:n_cells
        clear cell_data;
                
               
        for odor_n = 1:n_odors
            odor_ni = odor_list(odor_n);
            
            
            stim_fr_saved = [];
            stim_end_fr_saved = [];
            long_trace = [];
            prev_trace_length = 0;
            for odor_dur_n = 1:length(odor_dur_list);
                
                odor_dur_ni = odor_dur_list(odor_dur_n);                        %actual odor duration in s
                stim_end_fr = ceil(stim_frame + (odor_dur_ni./frame_time) );    %frame number when stim ends
                                
                tr_list = odor_trial_list_builder_20160317(stim_mat, prot_switch_trials, odor_ni, odor_dur_ni, block_n, an_trial_window, 1);
                                
                cell_data.traces(:, odor_n, odor_dur_n) = nanmean(squeeze(dff_data_mat(:, cell_n, tr_list, odor_ni) ), 2);
                cell_data.sds(:, odor_n, odor_dur_n) = nanstd(squeeze(dff_data_mat(:, cell_n, tr_list, odor_ni) ), [], 2);
                
                
                stim_fr_saved = [stim_fr_saved; stim_frame];
                stim_end_fr_saved = [stim_end_fr_saved; stim_end_fr];
                
                curr_trace = cell_data.traces(:, odor_n, odor_dur_n);
                extra_frames = 20;
                resp_trace = curr_trace((stim_frame-10):min([(stim_end_fr + extra_frames), length(curr_trace)]) );
                long_trace = [long_trace; resp_trace];
                
                stim_frs_saved(odor_dur_n, :) = [10, ((stim_end_fr - stim_frame) + 10)] + prev_trace_length;
                prev_trace_length = prev_trace_length + length(resp_trace);
                
            end
            
            cell_data.stim_start_frs = stim_fr_saved;
            cell_data.stim_end_frs = stim_end_fr_saved;
            
            
            %building large matrix with long ave traces for each cell in
            %this dataset
            big_matrix(:, cell_n, odor_n) = long_trace; 
            
            
        end
   
        %analysing current cell, curent odor
        %smoothing traces, moving window average
        traces = cell_data.traces;
        traces_s = filter([.2, .2, .2, .2, .2], 1, traces);
        
                
        cell_counter = cell_counter + 1;
        saved_cell_data{1, cell_counter} = cell_data;
    end
    %collating big matrices across datasets in this list
    bigger_matrixi = bigger_matrix;
    n_odors_orig = size(bigger_matrixi, 3); 
    n_cells_orig = size(bigger_matrixi, 2);      %n cells in all datasets upto current one
    n_frames_orig = size(bigger_matrixi, 1);     %n frames in all datasets upto current one
    n_odors_new = size(big_matrix, 3);
    n_cells_new = size(big_matrix, 2);           %n cells in current dataset
    n_frames_new = size(big_matrix, 1);          %n frames in current dataset
    bigger_matrix = zeros(max([n_frames, n_frames_new]), (n_cells_orig + n_cells_new), n_odors) + nan;
    if direc_counter > 1
        bigger_matrix(1:n_frames_orig, 1:n_cells_orig, 1:n_odors_new) = bigger_matrixi;
        bigger_matrix(1:n_frames_new, (n_cells_orig + 1):(n_cells_orig + n_cells_new), 1:n_odors_new) = big_matrix;
    elseif direc_counter == 1
        bigger_matrix(1:n_frames_new, (n_cells_orig + 1):(n_cells_orig + n_cells_new), 1:n_odors_new) = big_matrix;
    else
    end
    
end
fclose(fid);

%Clustering Analysis
%clustering cells as per responses to each odor separately

n_cells = size(bigger_matrix, 2);
for odor_n = 1:n_odors
    odor_ni = odor_list(odor_n);
    response_matrix = squeeze(bigger_matrix(:, :, odor_n));
    
    %------------------------------------
    X = response_matrix';       %activity data matrix - cells x frames, with trials concatenated

    thr = 0.8;              %threshold of kmeans++ runs a pair of cells must co-occur in to be considered for clustering
    no_iters = 500;         %no of times meta_k_means_tank calls k_means_plusplus; set to 100 for preliminary corrcut scan
    clustat = [];
    un_clust = [];
    
    ccorr = corrcoef(response_matrix, 'rows', 'pairwise');
   
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
    [a b c] = meta_k_means_tank1(X, 'correlation', corrcut, thr, no_iters);
    
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
    beep


    response_matrix_o = [];
    %accounting for un-clustered cells
    c_list = find(c_list_f == 0);
    un_clust = [un_clust; length(c_list), n_cells];
    for c_no = 1:length(c_list)
        c_noi = c_list(c_no);
        
    end

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
        curr_trace_mat = response_matrix(:, c_list)';       %matrix of activity data for each cell that belongs to this cluster

        curr_trace_mat = [cr_list, curr_trace_mat];     %concatenating corrcoefs to activity traces to sort both together
        curr_trace_mat = sortrows(curr_trace_mat, -1);      %sorting rows by first column ie corrcoefs
        response_matrix_o = [response_matrix_o; curr_trace_mat(:, 2:end)];
    end

    response_matrix_o = response_matrix_o';
    
    %calculating new corrcoeff mat with re-arranged cells
    curr_color = color_vec(odor_ni, :);
    figure(1)
    imagesc(response_matrix_o', [0, 4]);
    add_stim_shading(1, stim_frs_saved, 0.20, curr_color)
    colormap(greymap)
    
    c_o = corrcoef(response_matrix_o, 'rows', 'pairwise');
    figure(2)
    imagesc(c_o)
    colormap('jet')
    
    c_o = corrcoef(response_matrix_o, 'rows', 'pairwise');
    figure(3)
    imagesc(ccorr)
    colormap('jet')
    %------------------------------------
    n_clusts = max(cell_gp_vec(:, 1));
    for clust_n = 1:n_clusts
        curr_clust_list = find(cell_gp_vec(:, 1) == (clust_n));
        curr_traces = response_matrix(:, curr_clust_list);
        curr_traces = curr_traces./repmat(max(curr_traces), size(curr_traces, 1), 1);
        
        ave_resp_traces(:, 1) = nanmedian(curr_traces, 2); 
        
        figure(4)
        plot(ave_resp_traces(:, clust_n), 'LineWidth', 2)
        add_stim_shading(4, stim_frs_saved, 0.20, curr_color)
        
        xlabel('time')
        ylabel('Normalized dF/F')
        disp(['Cluster size ' num2str(length(curr_clust_list)./size(response_matrix, 2))])
        del = input('press enter');
    end
    
    
    %plotting trace for un-clustered cells
    clust_n = 0;
    curr_clust_list = find(cell_gp_vec(:, 1) == (clust_n));
    curr_traces = response_matrix(:, curr_clust_list);
    curr_traces = curr_traces./repmat(max(curr_traces), size(curr_traces, 1), 1);

    ave_resp_trace(:, 1) = nanmedian(curr_traces, 2); 

    figure(4)
    plot(ave_resp_trace(:, 1), 'LineWidth', 2)
    add_stim_shading(4, stim_frs_saved, 0.20, curr_color)
    title('Ave trace for un-clustered cells')
    
    xlabel('time')
    ylabel('Normalized dF/F')
    disp(['Cluster size ' num2str(length(curr_clust_list)./size(response_matrix, 2))])
    del = input('press enter');
    
end
