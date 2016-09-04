clear all
close all

list_direc = ['D:\Data\CSHL\dataset_list_dopamine_20150421.txt']; %expt datasets

blocks2compare = [1, 2];            %Specify the two blocks between which to compare KC pop responses
treatment_odor = 2;                 %The odor paired with something that is to be compared with other odors

an_trial_window = 5;                %Specify the window of trials to use from each block, for each odor (eg. first 4 trials)    
fid = fopen(list_direc);
direc_counter = 0;

plot_summary_an_only = 1;       %1 - only plots analyses across datasets, 0 - also plots for individual datasets
plot_PC_trajs = 1;              %if plot_summary_an_only is 0, this comes into effect
summary_statistic = 1;          % 0 - use mean; 1 - use median

sig_cell_block_mat = [];

all_c_dists = [];
mean_c_dists_all = [];
c_dist_lengths = [];
all_block1_dists = [];
color_vec = [[.75, .35, .25]; [.50, .70, .30]; [.55, .35, .75]; [.35, .55, .75]; [.75, .55, .35]; [.75, .35, .55]];
n_cells_saved = [];
mean_vec_first_pooled = [];
mean_vec_last_pooled = [];
se_vec_first_pooled = [];
se_vec_last_pooled = [];
odor1_corrs = [];
odor2_corrs = [];
all_c_dists_cosine = [];

%figure property initialisation variables
plot_height = 200;
plot_width = 200;
axis_font_size = 8;
tick_length = [0.02, 0.05];
line_widths = 2;
marker_sizes = 1;
marker_sizes_f = 4;
box_l_width = 0.5;

centroid_dists_saved = [];


%loop to go through all experiment datasets listed in list file
while 1
    direc_counter = direc_counter + 1;
    direc = fgetl(fid);
       
    if ischar(direc) ~= 1
        break
    else
    end
    
    direc = [direc '\'];
        
    %replacing C: with D:
    a = findstr(direc, 'C:');
    direc(a) = 'D';

    
    
    %loading extracted raw fluorescence data matrices written by
    %raw_dff_extractor
    dataset = load([direc 'expt.mat']);
    dataset = dataset.data;
    raw_data_mat = load([direc 'expt_raw_traces.mat']);
    raw_data_mat = raw_data_mat.raw_data_mat;
    
    
    %calculating dF/F traces and creating the sparse, 4-D, nan-filled
    %dff_data_mat 
    [dff_data_mat, stim_mat, prot_switch_trials] = cal_dff_traces_20150710(raw_data_mat, dataset, list_direc);
    clear raw_data_mat
    odor_list = unique(stim_mat(:, 1) );
    n_frames = size(dff_data_mat, 1);
    n_cells = size(dff_data_mat, 2);
    n_trials = size(dff_data_mat, 3);
    stim_time = dataset(1).stim.stimLatency;
    frame_time = dataset(1).info.framePeriod;
    stim_frame = floor(stim_time./frame_time);
    
   
       
    %identifying sig responses on a single trial basis, and then sig
    %responder cells in any individual block
    [resp_areas, sig_trace_mat, sig_cell_mat, sig_cell_block_mati] = cal_sig_responses_20151020(dataset, dff_data_mat, stim_mat, prot_switch_trials, list_direc, an_trial_window);
    
    
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
        for block_n = 1:n_prot_blocks
            od_blk_window_trs = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, block_n, an_trial_window, 1);
            
            n_bad_trs = length(intersect(od_blk_window_trs, bad_tr_list));            %number of bad trials in current block
            if n_bad_trs > curr_window./2
                disp('More than half the trials missing for one odor for one block.')
                bad_od_blk_counter = bad_od_blk_counter + 1;
            else
            end
            
        end
    end
    if bad_od_blk_counter > 1
        
        %continue
    else
    end
    
    

    sig_cell_block_mat = [sig_cell_block_mat; sig_cell_block_mati];
    
       
    %collecting block1 sig_responder cells' areas under curve in all trials 
    if direc_counter == 1
        for odor_na = 1:length(odor_list)
            saved_sig_areas{1, odor_na} = [];
            saved_insig_areas{1, odor_na} = [];
            saved_sig_traces = [];
        end

    else
    end

    for odor_n = 1:length(odor_list)
        odor_ni = odor_list(odor_n);    %odor_ni is actual odor number delivered, as opposed to position in list of odors
        odor_tr_list = find(stim_mat(:, 1) == odor_ni);
        non_odor_tr_list = 1:1:n_trials;
        del = intersect(non_odor_tr_list, odor_tr_list);
        non_odor_tr_list(del) = [];
        
        sig_cell_list = find(sig_cell_mat(:, odor_ni) == 1);
        insig_cell_list = find(sig_cell_mat(:, odor_ni) == 0);
                
        
        
        odor_resps = zeros(length(sig_cell_list), size(resp_areas, 2)) + nan;
        odor_resps_insig = zeros(length(insig_cell_list), size(resp_areas, 2)) + nan;
        odor_resps(:, odor_tr_list) = resp_areas(sig_cell_list, odor_tr_list);
        odor_resps_insig(:, odor_tr_list) = resp_areas(insig_cell_list, odor_tr_list);
        
               
        saved_sig_areas{1, odor_n} = [saved_sig_areas{1, odor_n}; odor_resps];
        saved_insig_areas{1, odor_n} = [saved_insig_areas{1, odor_n}; odor_resps_insig];
             
    end
    
    
    %PLOTTING
    odor_list = unique(stim_mat(:, 1));             %list of odor stimuli used in this dataset 
    n_odors = length(odor_list);       
    
    
    pre_tr_list = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, 1, an_trial_window, 0);
    post_tr_list = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, length(prot_switch_trials), an_trial_window, 1);
    
    if plot_summary_an_only == 0
        %plotting trial-averaged response traces
        figure(1)
        for odor_ni = 1:n_odors
            odor_n = odor_list(odor_ni);
            odor_name = dataset(1).stim.odourNames(odor_n).odour;
            try
                subplot(ceil(sqrt(n_odors)), ceil(sqrt(n_odors)), odor_ni)
            catch
                keyboard
            end
            imagesc(squeeze(nanmean(dff_data_mat(:, :, :, odor_n), 3))', [0, 3]);
            title(['Trial-averaged traces for ' odor_name])
            xlabel('frame number')
            ylabel('cell number')

        end
        set(gcf, 'Color', 'w')
        a = colormap('bone');
        a = flipud(a);
        colormap(a)

        %Plotting response areas over all trials
        figure(2)
        for odor_ni = 1:n_odors
            odor_n = odor_list(odor_ni);
            odor_trs = find(stim_mat(:, 1) == odor_n);

            odor_name = dataset(1).stim.odourNames(odor_n).odour;
            subplot(ceil(sqrt(n_odors)), ceil(sqrt(n_odors)), odor_ni)
            imagesc(resp_areas(:, odor_trs), [0, 30])

            title(['Area under the peak for ' odor_name])
            xlabel('trial number')
            ylabel('cell number')
    
        end
        
        a = colormap('bone');
        a = flipud(a);
        colormap(a)
        
        
        %plotting response areas across only analysed trials, only for significantly responsive cells 
        figure(9)
        %tr_list = [3:10, 13:20];
        for odor_ni = 1:n_odors
            odor_n = odor_list(odor_ni);
            odor_trs = find(stim_mat(:, 1) == odor_n);
            del = find(odor_trs > 20);
            odor_trs(del) = [];
            sig_cell_list_temp = find(sig_cell_mat(:, odor_n) == 1);

            odor_name = dataset(1).stim.odourNames(odor_n).odour;
            subplot(ceil(sqrt(n_odors)), ceil(sqrt(n_odors)), odor_ni)
            imagesc(resp_areas(sig_cell_list_temp, odor_trs), [0, 30])

            title(['Area under the peak for ' odor_name])
            xlabel('trial number')
            ylabel('cell number')
    
        end
        
        
        
        a = colormap('bone');
        a = flipud(a);
        colormap(a)
             
        
        set(gcf, 'Color', 'w')
    
        %Plotting mean response traces over trials, separated by protocol block
        fig_handle3 = figure(3)
        subplot_counter = 0;
        for prot_block_n = 1:n_prot_blocks
            if prot_block_n == 1
                tr_start = 1;                      %first trial of current protocol block
            else
                tr_start = prot_switch_trials(prot_block_n - 1);
            end

            tr_end = prot_switch_trials(prot_block_n) - 1;  %last trial of current protocol block
            
            a = colormap('bone');
            a = flipud(a);
            colormap(a)
            
            for odor_ni = 1:n_odors
                odor_n = odor_list(odor_ni);
                subplot_counter = subplot_counter + 1;
                fig_handle = figure(subplot_counter + 10);
                odor_name = dataset(1).stim.odourNames(odor_n).odour;
                %subplot(n_prot_blocks, n_odors, subplot_counter)
                
%                 if prot_block_n == 1
%                     imagesc(squeeze(nanmean(dff_data_mat(:, :, (tr_end-8):tr_end, odor_n), 3))', [0, 3]);
%                 elseif prot_block_n == 3
%                     imagesc(squeeze(nanmean(dff_data_mat(:, :, tr_start:(tr_start + 8), odor_n ), 3))', [0, 3]);
%                 else
                    imagesc(squeeze(nanmean(dff_data_mat(:, :, tr_start:tr_end, odor_n), 3))', [0, 3]);
%                 end
                    
                title(['Trial Block ' int2str(prot_block_n) '; ' odor_name])
                xlabel('frame number')
                ylabel('cell number')
                set(fig_handle, 'Position', [100, 100, plot_height, plot_width]);
                set(gca,'FontSize', axis_font_size, 'ticklength', tick_length, 'linewidth', box_l_width)
                
                a = colormap('bone');
                a = flipud(a);
                colormap(a)
            end
        end
        set(gcf, 'Color', 'w')
        
        
        if n_cells < 3
            continue
        else
            
            if length(good_tr_list) < 3
                continue
            else
            end
            
            %PCA for plotting pop response trajectories
            nPCs = 2;
            [dff_all_PCs, areas_all_PCs ] = PCA_calculator(resp_areas, nPCs, dff_data_mat, odor_list);
            
            area_mat_for_plotting = zeros(1, nPCs, n_trials) + nan;
            area_mat_for_plotting(1, :, :) = areas_all_PCs;
            
            if plot_PC_trajs == 1 
                %plotting response time-courses in PC space
                for plot_odor_n = 1:n_odors
                    odor_ni = odor_list(plot_odor_n);
                    pre_tr_list = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, 1, an_trial_window, 0);
                    post_tr_list = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, length(prot_switch_trials), an_trial_window, 1);

                    plot_data = (dff_all_PCs((stim_frame - 2):(stim_frame + 8), :, pre_tr_list, plot_odor_n));
                    trajectory_plotter(plot_data, 4, [], 0, (color_vec(plot_odor_n, :) + 0.15), 0, 0, marker_sizes_f);

                    del = input('press enter for post pairing trials');
                    plot_data = (dff_all_PCs((stim_frame - 2):(stim_frame + 8), :, post_tr_list, plot_odor_n));
                    figure(4)
                    trajectory_plotter(plot_data, 4, [], 0, (color_vec(plot_odor_n, :) - 0.2), 0, 0, marker_sizes_f);


                end
                set(gcf, 'Position', [100, 100, plot_height, plot_width]);
                set(gca,'FontSize', axis_font_size, 'ticklength', tick_length, 'linewidth', box_l_width)

            
                %plotting response areas in PC space
                for plot_odor_n = 1:n_odors
                    odor_ni = odor_list(plot_odor_n);
                    pre_tr_list = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, 1, an_trial_window, 0);
                    post_tr_list = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, length(prot_switch_trials), an_trial_window, 1);

                    plot_data = (area_mat_for_plotting(1, :, pre_tr_list));
                    fig_handle4 = figure(4);
                    trajectory_plotter(plot_data, 5, [], 1, (color_vec(plot_odor_n, :) + 0.15), 0, 0, marker_sizes_f);


                    del = input('press enter for post pairing trials');
                    plot_data = (area_mat_for_plotting(1, :, post_tr_list));
                    fig_handle5 = figure(5);
                    trajectory_plotter(plot_data, 5, [], 1, (color_vec(plot_odor_n, :) - 0.2), 0, 0, marker_sizes_f);


                end
                set(fig_handle4, 'Position', [100, 100, plot_height, plot_width]);
                set(gca,'FontSize', axis_font_size, 'ticklength', tick_length, 'linewidth', box_l_width)
                xlabel('PC 1')
                ylabel('PC 2')
                set(fig_handle5, 'Position', [100, 100, plot_height, plot_width]);
                set(gca,'FontSize', axis_font_size, 'ticklength', tick_length, 'linewidth', box_l_width)
                xlabel('PC 1')
                ylabel('PC 2')
            else
            end
            
        
        end

    else
    end  


    
    
    %calculating population representation distances for current fly
    if n_cells > 3
        
        for odor_n = 1:length(odor_list)
            odor_ni = odor_list(odor_n);
            %curr_sig_cells = find(sig_cell_mat(:, odor_ni) == 1);
            curr_sig_cells = 1:n_cells;                     %changed our minds and wanted to include all cells
            
            an_tr_list1 = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, 1, an_trial_window, 0);        %list of trial_window number of trials for this odor from the end of block1 
            an_tr_list2 = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, length(prot_switch_trials), an_trial_window, 1);        %list of trial_window number of trials for this odor from the end of the last block
            
            block1_areas = resp_areas(curr_sig_cells, an_tr_list1);
            c_dist_lengths = [c_dist_lengths; length(curr_sig_cells)];
            block1_centroid = nanmean(block1_areas, 2);
            block1_coords = [block1_centroid, block1_areas];
            block1_dists = squareform(pdist(block1_coords'));
            block1_dists_cosine = squareform(pdist(block1_coords', 'cosine'));
            
            block1_dists = (block1_dists(1, 2:end));                    %distances of block1 trials from block1 centroid
            block1_dists_cosine = (block1_dists_cosine(1, 2:end));                    %distances of block1 trials from block1 centroid
            
            block1_centroid_dist = squareform(pdist([block1_centroid, zeros(size(block1_centroid, 1), 1)]'));
            block1_centroid_dist_cosine = squareform(pdist([block1_centroid, zeros(size(block1_centroid, 1), 1)]', 'cosine'));
            block1_centroid_dist = block1_centroid_dist(1, 2);          %distance of block1 centroid from origin; used later for normalisation
                        
            block1_dists = block1_dists./block1_centroid_dist;          %normalising distance to pre-pairing centroid
            block1_dists_cosine = block1_dists_cosine;          %not normalised as cosine dist from origin is not defined hence centroid dist not defined
            block1_dists_mat(:, odor_n) = block1_dists;
            block1_dists_mat_cosine(:, odor_n) = block1_dists_cosine;
            
            
            block2_areas = resp_areas(curr_sig_cells, an_tr_list2);
            block2_centroid = nanmean(block2_areas, 2);
            block2_centroid_dist = squareform(pdist([block2_centroid, zeros(size(block2_centroid, 1), 1)]'));
            block2_centroid_dist = block2_centroid_dist(1, 2);          %distance of block2 centroid from origin; 
            centroid_dists_saved = [centroid_dists_saved; [block1_centroid_dist, block2_centroid_dist]];        %saving distances of all block1 centroids from origin
            
            coords_mat = [block1_centroid, block2_areas];
            c_dists = squareform(pdist(coords_mat'));                   %symmetric distance matrix of each trial from the centroid
            c_dists_cosine = squareform(pdist(coords_mat', 'cosine'));                   %symmetric distance matrix of each trial from the centroid
            
            c_dists = c_dists(1, 2:end)./block1_centroid_dist;          %list of distances of each block2 trial from the block1 centroid
            c_dists_cosine = c_dists_cosine(1, 2:end);          %not normalising as cosine dist from origin is not defined, hence centroid dist from origin is not defined
            
            c_dist_mat(:, odor_n) = c_dists;
            c_dist_mat_cosine(:, odor_n) = c_dists_cosine;
            
            mean_c_dists(1, odor_n) = nanmean(c_dists);
            
            
            %calculating population representation correlations for current fly.
            X_block2_areas = [block1_centroid, block2_areas];     %concatenating pre-training centroid to post-training trial responses
            corrs = corrcoef(X_block2_areas);
            corrs = corrs(1:5, 1);

            if odor_n == 1
                odor1_corrs = [odor1_corrs; corrs];
            elseif odor_n == 2
                odor2_corrs = [odor2_corrs; corrs];
            else
            end
            
            
            
        end
        all_c_dists = [all_c_dists; c_dist_mat];
        all_c_dists_cosine = [all_c_dists_cosine; c_dist_mat_cosine];
        clear c_dist_mat
        clear c_dist_mat_cosine
        
        mean_c_dists_all = [mean_c_dists_all; mean_c_dists];
        clear mean_c_dists
        
        all_block1_dists = [all_block1_dists; block1_dists_mat];
        clear block1_dists_mat
        
    
    
    else
    end
    
    %keyboard
    
    figure(4)
    close gcf

    figure(5)
    close gcf
    
end
fclose(fid);





%PLOTTING/ANALYSING across datasets
sorted_mean_vecs_first = [];
sorted_mean_vecs_last = [];
for odor_n = 1:length(odor_list)
    odor_ni = odor_list(odor_n);
    odor_name = dataset(1).stim.odourNames(odor_n).odour;
    odor_resp_mat = saved_sig_areas{1, odor_n};
    odor_resp_mat_insig = saved_insig_areas{1, odor_n};
    curr_plot_color = color_vec(odor_n, :);
        
    %Reality check: does average activity of KC population change across
    %protocol blocks?
    
    fig_handle = figure(1);
    if summary_statistic == 0
        mean_vec = nanmean(odor_resp_mat(:, : ) );                               %discarding trial 1 for each odor
        se_vec = nanstd(odor_resp_mat(:, :))./sqrt(size(odor_resp_mat, 1));      %discarding trial 1
        
    elseif summary_statistic == 1
        mean_vec = nanmedian(odor_resp_mat(:, :) );                             %discarding trial 1
        se_vec = mad(odor_resp_mat(:, :))./sqrt(size(odor_resp_mat, 1));      %discarding trial 1
    else
    end
    
    %resizing figure
    
        
    errorbar([1:length(mean_vec)], mean_vec, se_vec, 'Color', curr_plot_color, 'LineWidth', 1)
    hold on
    plot(mean_vec, ['O' curr_plot_color], 'markerfacecolor', curr_plot_color, 'markeredgecolor', curr_plot_color)
    
    set(fig_handle, 'Position', [100, 100, plot_height, plot_width]);
    set(gca,'FontSize', axis_font_size, 'ticklength', tick_length, 'linewidth', box_l_width)
    
    xlabel('Trial number')
    if summary_statistic == 0
        ylabel('Mean response area across cells')
    elseif summary_statistic == 1
        ylabel('Median response area across cells')
    else
    end
    
    
    
    %plotting mean and error bar for each cell's responses in first v/s
    %second of selected blocks
    fig_handle2 = figure(2);
    an_tr_list1 = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, blocks2compare(1), an_trial_window, 0);        %list of trial_window number of trials for this odor from the end of block1 
    an_tr_list2 = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, blocks2compare(2), an_trial_window, 1);        %list of trial_window number of trials for this odor from the end of the last block
    an_tr_list = [an_tr_list1; an_tr_list2];
    
    if summary_statistic == 0
        mean_vec_first = nanmean(odor_resp_mat(:, an_tr_list1), 2);
        mean_vec_last = nanmean(odor_resp_mat(:, an_tr_list2), 2);
        se_vec_first = nanstd(odor_resp_mat(:, an_tr_list1), [], 2)./sqrt(length(an_tr_list));    
        se_vec_last = nanstd(odor_resp_mat(:, an_tr_list2), [], 2)./sqrt(length(an_tr_list));    
        
        %calculating separately for insignificantly responsive cells
        mean_vec_first_insig = nanmean(odor_resp_mat_insig(:, an_tr_list1), 2);
        mean_vec_last_insig = nanmean(odor_resp_mat_insig(:, an_tr_list2), 2);
        se_vec_first_insig = nanstd(odor_resp_mat_insig(:, an_tr_list1), [], 2)./sqrt(length(an_tr_list));    
        se_vec_last_insig = nanstd(odor_resp_mat_insig(:, an_tr_list2), [], 2)./sqrt(length(an_tr_list));    
        
    elseif summary_statistic == 1
        mean_vec_first = nanmedian(odor_resp_mat(:, an_tr_list1), 2);
        mean_vec_last = nanmedian(odor_resp_mat(:, an_tr_list2), 2);
        se_vec_first = mad(odor_resp_mat(:, an_tr_list1), [], 2)./sqrt(length(an_tr_list));    
        se_vec_last = mad(odor_resp_mat(:, an_tr_list2), [], 2)./sqrt(length(an_tr_list));    
        
        %calculating separately for insignificantly responsive cells
        mean_vec_first_insig = nanmedian(odor_resp_mat_insig(:, an_tr_list1), 2);
        mean_vec_last_insig = nanmedian(odor_resp_mat_insig(:, an_tr_list2), 2);
        se_vec_first_insig = mad(odor_resp_mat_insig(:, an_tr_list1), [], 2)./sqrt(length(an_tr_list));    
        se_vec_last_insig = mad(odor_resp_mat_insig(:, an_tr_list2), [], 2)./sqrt(length(an_tr_list));    
        
        
    else
    end
        
    %dropping padding nans
    del = isnan(mean_vec_first);
    mean_vec_first(del) = [];
    se_vec_first(del) = [];
    
    mean_vec_last(del) = [];
    se_vec_last(del) = [];
    
    
    %cumulative hist of response sizes pre treatment and post treatment for
    %each odor
    fig_handle4 = figure(4);
    bin_vec = -0.6:0.3:3;
    hist_first = hist(mean_vec_first, bin_vec);
    hist_first = hist_first./sum(hist_first);
    hist_first = cumsum(hist_first);
    
    hist_last = hist(mean_vec_last, bin_vec);
    hist_last = hist_last./sum(hist_last);
    hist_last = cumsum(hist_last);
    
    plot(bin_vec, hist_first, 'Color', (curr_plot_color - 0.2), 'LineWidth', 1)
    hold on
    plot(bin_vec, hist_last, 'Color', (curr_plot_color + 0.2), 'LineWidth', 1)
    
    set(fig_handle4, 'Position', [100, 100, plot_height, plot_width]);
    set(gca,'FontSize', axis_font_size, 'ticklength', tick_length, 'linewidth', box_l_width)
    
    xlabel('Response size');
    ylabel('Cumulative fraction of cells')
    
    fig_handle2 = figure(2);
    
    %making two plots, one without non-significant cell resposnes and
    %another with them included
    for copy_n = 1:4
        if copy_n == 1
            fig_handle2 = figure(2)
        elseif copy_n == 2
            fig_handle2 = figure(10)
        else
        end
        
        errorbar(mean_vec_first, mean_vec_last, se_vec_last, 'O', 'MarkerFaceColor', curr_plot_color, 'MarkerEdgeColor', curr_plot_color, 'MarkerSize', marker_sizes, 'Color', curr_plot_color)
        hold on
        h = herrorbar(mean_vec_first, mean_vec_last, se_vec_first, ['.' curr_plot_color])
        set(h, 'Color', curr_plot_color)
        xlabel('Single cell responses in first block')
        ylabel('Single cell responses in last block')
        %set(gcf, 'Color ', 'w')

        %fitting lines to responses in block1 v/s last block
        line_pars = polyfit(mean_vec_first, mean_vec_last, 1);
        linex = [0:1:max([mean_vec_first; mean_vec_last])];
        liney = polyval(line_pars, linex);


        plot(linex, liney, '-', 'Color', curr_plot_color, 'LineWidth', line_widths./2)

        
        %for second copy of figure with responses from insignificantly
        %responsive cells also included
        if copy_n == 2
            curr_plot_color = curr_plot_color + 0.2;
            errorbar(mean_vec_first_insig, mean_vec_last_insig, se_vec_last_insig, 'O', 'MarkerFaceColor', curr_plot_color, 'MarkerEdgeColor', curr_plot_color, 'MarkerSize', marker_sizes, 'Color', curr_plot_color)
            hold on
            h = herrorbar(mean_vec_first_insig, mean_vec_last_insig, se_vec_first_insig, ['.' curr_plot_color])
            set(h, 'Color', curr_plot_color)
            xlabel('Single cell responses in first block')
            ylabel('Single cell responses in last block')
            %set(gcf, 'Color ', 'w')

            %fitting lines to responses in block1 v/s last block
            line_pars = polyfit(mean_vec_first_insig, mean_vec_last_insig, 1);
            linex = [0:1:max([mean_vec_first_insig; mean_vec_last_insig])];
            liney = polyval(line_pars, linex);


            plot(linex, liney, '-', 'Color', curr_plot_color, 'LineWidth', line_widths./2)
            curr_plot_color = curr_plot_color - 0.2;

        else
        end
        
        
        
        %plotting slope - 1 reference line for figure 2
        axisi = axis;
        del_max = max(axisi);
        del_min = min(axisi);
        axis([del_min, del_max, del_min, del_max])
        plot([del_min, del_max], [del_min, del_max], '--', 'Color', [0.60, 0.60, 0.60])
        set(fig_handle2, 'Position', [100, 100, plot_height, plot_width]);
        set(gca,'FontSize', axis_font_size, 'ticklength', tick_length, 'linewidth', box_l_width)
       
        
    end

    
    
    % checking 
      
        
    %STATISTICS (single-cell level analysis)
    %checking if response dists are normal
    %tabulating counts of significant cells
    T = [];

    odor_ni = odor_list(odor_n);
    odor_name = dataset(1).stim.odourNames(odor_ni).odour;
    
    odor_num(odor_n) = odor_n;
    block1_sig_cells(odor_n, 1) = nansum(sig_cell_block_mat(:, odor_ni, blocks2compare(1) ));
    block2_sig_cells(odor_n, 1) = nansum(sig_cell_block_mat(:, odor_ni, blocks2compare(2) ));
    del = sig_cell_block_mat(:, :, blocks2compare(2)) - sig_cell_block_mat(:, :, blocks2compare(1));
    block1not2_sig_cells(odor_n, 1) = length(find(del(:, odor_ni) == -1) );
    block2not1_sig_cells(odor_n, 1) = length(find(del(:, odor_ni) == 1) );
    block1n2_sig_cells(odor_n, 1) = block1_sig_cells(odor_n) - block1not2_sig_cells(odor_n);
    T_line = [{odor_name}, block1_sig_cells(odor_n, 1), block2_sig_cells(odor_n, 1), block1not2_sig_cells(odor_n, 1), block2not1_sig_cells(odor_n, 1), block1n2_sig_cells(odor_n, 1)];
    T = [T; T_line];
    T = [{'Odor'}, {'Block1'}, {'Block2'}, {'Block1not2'}, {'Block2not1'}, {'Block1and2'}; T];
    disp('summary significant cell counts')
    disp(T)

    data_to_save{1, odor_n} = [mean_vec_first, mean_vec_last];
    data_to_save_ses{1, odor_n} = [se_vec_first, se_vec_last];
    
    
    %Plotting 
    curr_odor_data = [mean_vec_first, mean_vec_last];
    odor_diff = curr_odor_data(:, 2) - curr_odor_data(:, 1);   %calculating change in dF/F response sizes across compared prot_blocks
    fig_handle3 = figure(3);
    n_bins = 12;
    bin_vec = -1:(4./n_bins):3;
    [odor_hist(:, 1), hist_centers(:, 1)] = hist(odor_diff, bin_vec);
    odor_hist(:, 1) = (odor_hist(:, 1)./sum(odor_hist(:, 1)) ).*100;           %normalising numbers of cells in each bin
    plot(hist_centers(:, 1), odor_hist(:, 1), 'Color', curr_plot_color, 'LineWidth', line_widths)
    hold on
    xlabel('Change in response size');
    ylabel('Percentage of cells');
    set(fig_handle3, 'Position', [100, 100, plot_height, plot_width]);
    set(gca,'FontSize', axis_font_size, 'ticklength', tick_length, 'linewidth', box_l_width)
  
    
    %looking at cells with large increases in response size
    big_changei = find(odor_diff > 1);       %looking at cells whose responses have more than doubled
    small_changei = find(odor_diff <= 1);
    big_ses = se_vec_last(big_changei);
    small_ses = se_vec_last(small_changei);
    
    mean_vecs = [mean_vec_first, mean_vec_last];
    mean_vecs_sr = sortrows(mean_vecs);
    
    fig_handle15 = figure(15)
    bar(mean_vecs_sr(:, 1))

    set(fig_handle15, 'Position', [100, 100, plot_height, plot_width]);
    set(gca,'FontSize', axis_font_size, 'ticklength', tick_length, 'linewidth', box_l_width)
    xlabel('Cell number')
    ylabel('Mean response size in Block1')
    
    
    fig_handle16 = figure(16)
    bar(mean_vecs_sr(:, 2))
    
    set(fig_handle16, 'Position', [100, 100, plot_height, plot_width]);
    set(gca,'FontSize', axis_font_size, 'ticklength', tick_length, 'linewidth', box_l_width)
    xlabel('cell number')
    ylabel('Mean response size in Block2')

    keyboard
end



%plotting N-dimensional distances
fig_handle5 = figure(5);
%boxplot(mean_c_dists_all, 'boxstyle', 'filled', 'colors', color_vec)
boxplot(all_c_dists, 'boxstyle', 'filled', 'colors', color_vec, 'widths', 0.2, 'symbol', 'o', 'medianstyle', 'target')
%plotting individual points with random jitter in x
hold on

for odor_n = 1:length(odor_list)
    curr_color = color_vec(odor_n, :);
    rand_vec = (rand(1, size(all_c_dists, 1)).*0.5 + odor_n - 0.25);
    
    plot(rand_vec, all_c_dists, 'o', 'Color', curr_color, 'MarkerSize', marker_sizes_f)
    
end

set(fig_handle5, 'Position', [100, 100, plot_height, plot_width]);
set(gca,'FontSize', axis_font_size, 'ticklength', tick_length, 'linewidth', box_l_width)
ylabel('Distances from pre-training centroid')
xlabel('Filler axis label')

%distances of post-pairing trials from pre-pairing centroid
disp('distances of post-pairing trials from pre-pairing centroid')
[h, p] = ttest2(all_c_dists(:, 1), all_c_dists(:, 2))

disp('cosine distances of post-pairing trials from pre-pairing centroid')
[h, p] = ttest2(all_c_dists_cosine(:, 1), all_c_dists_cosine(:, 2))


%distances of pre-pairing centroids from origin
%Can't really use this since each distance is measured in a different
%dimensionality. Not normalised because these distances are what are used
%for the normalisation in the first place.
[h, p] = ttest2(centroid_dists_saved(1:(length(centroid_dists_saved)./2), 1), centroid_dists_saved( ((length(centroid_dists_saved)./2)+1):length(centroid_dists_saved), 1));

%distances of post-pairing centroids from origin
centroid_dists_saved(:, 2) = centroid_dists_saved(:, 2)./centroid_dists_saved(:, 1);        %normalising
disp('distances of post-pairing centroids from origin')
[h, p] = ttest2(centroid_dists_saved(1:(length(centroid_dists_saved)./2), 2), centroid_dists_saved( ((length(centroid_dists_saved)./2)+1):length(centroid_dists_saved), 2))



%correlation measure of distance of post-pairing trials from pre-centroid
disp('correlations of post-pairing trials from pre-pairing centroids')
[h, p] = ttest2(odor1_corrs, odor2_corrs)





