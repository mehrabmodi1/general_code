clear all
close all
%list_direc = ['C:\Data\CSHL\dataset_list_Toshi_KC_DN_Led_20150708.txt']; %expt datasets
list_direc = ['C:\Data\CSHL\dataset_KC_DN_pairing_control_20150828.txt']; %control datasets  

fid = fopen(list_direc);
direc_counter = 0;

plot_summary_an_only = 1;       %1 - only plots analyses across datasets, 0 - also plots for individual datasets
plot_PC_trajs = 0;              %if plot_summary_an_only is 0, this comes into effect
blk_trial_window = 4;               %uses the last 'trial_window' trials in block1 and the first 'trial_window' trials in blockLast for analysis, plotting. If trial_window is empty, the entire dataset is analysed.
summary_statistic = 1;          % 0 - use mean; 1 - use median

sig_cell_block_mat = [];
all_c_dists = [];
mean_c_dists_all = [];
c_dist_lengths = [];
all_block1_dists = [];
color_vec = [[.75, .35, .25]; [.50, .70, .30]; [.55, .35, .75]; [.35, .55, .75]; [.75, .55, .35]; [.75, .35, .55]];
n_cells_saved = [];


%figure property initialisation variables
plot_height = 300;
plot_width = 300;
axis_font_size = 15;
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
    
    n_frames = size(dff_data_mat, 1);
    n_cells = size(dff_data_mat, 2);
    n_trials = size(dff_data_mat, 3);
    stim_time = dataset(1).stim.stimLatency;
    frame_time = dataset(1).info.framePeriod;
    stim_frame = floor(stim_time./frame_time);
    
    n_cells_saved = [n_cells_saved; n_cells];
       
    %identifying sig responses on a single trial basis, and then sig
    %responder cells in any individual block
    [resp_areas, sig_trace_mat, sig_cell_mat, sig_cell_block_mati] = cal_sig_responses_20150825(dataset, dff_data_mat, stim_mat, prot_switch_trials, list_direc);
    sig_cell_block_mat = [sig_cell_block_mat; sig_cell_block_mati];
    
       
    %collecting block1 sig_responder cells' areas under curve in all trials 
    odor_list = unique(stim_mat(:, 1) );
    for odor_n = 1:length(odor_list)
        odor_ni = odor_list(odor_n);    %odor_ni is actual odor number delivered, as opposed to position in list of odors
        odor_tr_list = find(stim_mat(:, 1) == odor_ni);
        non_odor_tr_list = 1:1:n_trials;
        del = intersect(non_odor_tr_list, odor_tr_list);
        non_odor_tr_list(del) = [];
        
        sig_cell_list = find(sig_cell_mat(:, odor_ni) == 1);
        
        if direc_counter == 1
            for odor_na = 1:length(odor_list)
                saved_sig_areas{1, odor_na} = [];
                saved_sig_traces = [];
            end
            
        else
        end
        
        odor_resps = zeros(length(sig_cell_list), size(resp_areas, 2)) + nan;
        odor_resps(:, odor_tr_list) = resp_areas(sig_cell_list, odor_tr_list);
        
               
        saved_sig_areas{1, odor_n} = [saved_sig_areas{1, odor_n}; odor_resps];
                
    end
    
    
    %PLOTTING
    odor_list = unique(stim_mat(:, 1));             %list of odor stimuli used in this dataset 
    n_odors = length(odor_list);       
    n_prot_blocks = length(prot_switch_trials);     %number of protocol blocks in this dataset
    
    pre_tr_list = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, 1, blk_trial_window, 0);
    post_tr_list = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, length(prot_switch_trials), blk_trial_window, 1);
    
    if plot_summary_an_only == 0
        %plotting trial-averaged response traces
        figure(1)
        for odor_ni = 1:n_odors
            odor_n = odor_list(odor_ni);
            odor_name = dataset(1).stim.odourNames(odor_n).odour;
            subplot(ceil(sqrt(n_odors)), ceil(sqrt(n_odors)), odor_ni)
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
        tr_list = [3:10, 13:20];
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
        figure(3)
        subplot_counter = 0;
        for prot_block_n = 1:n_prot_blocks
            if prot_block_n == 1
                tr_start = 1;                      %first trial of current protocol block
            else
                tr_start = prot_switch_trials(prot_block_n - 1) + 1;
            end

            tr_end = prot_switch_trials(prot_block_n);  %last trial of current protocol block
            
            a = colormap('bone');
            a = flipud(a);
            colormap(a)
            
            for odor_ni = 1:n_odors
                odor_n = odor_list(odor_ni);
                subplot_counter = subplot_counter + 1;

                odor_name = dataset(1).stim.odourNames(odor_n).odour;
                subplot(n_prot_blocks, n_odors, subplot_counter)
                
                if prot_block_n == 1
                    imagesc(squeeze(nanmean(dff_data_mat(:, :, (tr_end-8):tr_end, odor_n), 3))', [0, 3]);
                elseif prot_block_n == 3
                    imagesc(squeeze(nanmean(dff_data_mat(:, :, tr_start:(tr_start + 8), odor_n ), 3))', [0, 3]);
                else
                    imagesc(squeeze(nanmean(dff_data_mat(:, :, tr_start:tr_end, odor_n), 3))', [0, 3]);
                end
                    
                title(['Trial Block ' int2str(prot_block_n) '; ' odor_name])
                xlabel('frame number')
                ylabel('cell number')

            end
        end
        set(gcf, 'Color', 'w')
        
        

        if n_cells < 3
            continue
        else
            %PCA for plotting pop response trajectories
            nPCs = 2;
            [dff_all_PCs, areas_all_PCs ] = PCA_calculator(resp_areas, nPCs, dff_data_mat, odor_list, blk_trial_window);
            
            area_mat_for_plotting = zeros(1, nPCs, n_trials) + nan;
            area_mat_for_plotting(1, :, :) = areas_all_PCs;
            
            if plot_PC_trajs == 1 
                %plotting response time-courses in PC space
                for plot_odor_n = 1:n_odors
                    odor_ni = odor_list(plot_odor_n);
                    pre_tr_list = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, 1, blk_trial_window, 0);
                    post_tr_list = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, length(prot_switch_trials), blk_trial_window, 1);

                    plot_data = (dff_all_PCs((stim_frame - 2):(stim_frame + 8), :, pre_tr_list, plot_odor_n));
                    trajectory_plotter(plot_data, 4, [], 0, (color_vec(plot_odor_n, :) + 0.15), 0, 0);

                    del = input('press enter for post pairing trials');
                    plot_data = (dff_all_PCs((stim_frame - 2):(stim_frame + 8), :, post_tr_list, plot_odor_n));
                    figure(4)
                    trajectory_plotter(plot_data, 4, [], 0, (color_vec(plot_odor_n, :) - 0.2), 0, 0);


                end
                set(gcf, 'Position', [100, 100, plot_height, plot_width]);
                set(gca,'FontSize', axis_font_size)

            
                %plotting response areas in PC space
                for plot_odor_n = 1:n_odors
                    odor_ni = odor_list(plot_odor_n);
                    pre_tr_list = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, 1, blk_trial_window, 0);
                    post_tr_list = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, length(prot_switch_trials), blk_trial_window, 1);

                    plot_data = (area_mat_for_plotting(1, :, pre_tr_list));
                    fig_handle4 = figure(4);
                    trajectory_plotter(plot_data, 5, [], 2, (color_vec(plot_odor_n, :) + 0.15), 0, 0);


                    del = input('press enter for post pairing trials');
                    plot_data = (area_mat_for_plotting(1, :, post_tr_list));
                    fig_handle5 = figure(5);
                    trajectory_plotter(plot_data, 5, [], 2, (color_vec(plot_odor_n, :) - 0.2), 0, 0);


                end
                set(fig_handle4, 'Position', [100, 100, plot_height, plot_width]);
                set(gca,'FontSize', axis_font_size)
                xlabel('PC 1')
                ylabel('PC 2')
                set(fig_handle5, 'Position', [100, 100, plot_height, plot_width]);
                set(gca,'FontSize', axis_font_size)
                xlabel('PC 1')
                ylabel('PC 2')
            else
            end
            
        
        end

      keyboard
      figure(4)
      close gcf
      
      figure(5)
      close gcf
    else
    end  


    
    
    %calculating population representation distances for current fly
    if n_cells > 3
        
        for odor_n = 1:length(odor_list)
            odor_ni = odor_list(odor_n);
            curr_sig_cells = find(sig_cell_mat(:, odor_ni) == 1);
            
            an_tr_list1 = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, 1, blk_trial_window, 0);        %list of trial_window number of trials for this odor from the end of block1 
            an_tr_list2 = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, length(prot_switch_trials), blk_trial_window, 1);        %list of trial_window number of trials for this odor from the end of the last block
            
            block1_areas = resp_areas(curr_sig_cells, an_tr_list1);
            c_dist_lengths = [c_dist_lengths; length(curr_sig_cells)];
            block1_centroid = nanmean(block1_areas, 2);
            block1_coords = [block1_centroid, block1_areas];
            block1_dists = squareform(pdist(block1_coords'));
            block1_dists = (block1_dists(1, 2:end));                    %distances of block1 trials from block1 centroid
            
            block1_centroid_dist = squareform(pdist([block1_centroid, zeros(size(block1_centroid, 1), 1)]'));
            block1_centroid_dist = block1_centroid_dist(1, 2);          %distance of block1 centroid from origin; used later for normalisation
            
            block1_dists = block1_dists./block1_centroid_dist;          %normalising
            block1_dists_mat(:, odor_n) = block1_dists;
            
            
            block2_areas = resp_areas(curr_sig_cells, an_tr_list2);
            block2_centroid = nanmean(block2_areas, 2);
            block2_centroid_dist = squareform(pdist([block2_centroid, zeros(size(block2_centroid, 1), 1)]'));
            block2_centroid_dist = block2_centroid_dist(1, 2);          %distance of block2 centroid from origin; 
            centroid_dists_saved = [centroid_dists_saved; [block1_centroid_dist, block2_centroid_dist]];        %saving distances of all block1 centroids from origin
            
            coords_mat = [block1_centroid, block2_areas];
            c_dists = squareform(pdist(coords_mat'));                   %symmetric distance matrix of each trial from the centroid
            
            c_dists = c_dists(1, 2:end)./block1_centroid_dist;          %list of distances of each block2 trial from the block1 centroid
            
            c_dist_mat(:, odor_n) = c_dists;  
            
            mean_c_dists(1, odor_n) = nanmean(c_dists);
            
        end
        all_c_dists = [all_c_dists; c_dist_mat];
        clear c_dist_mat
        
        mean_c_dists_all = [mean_c_dists_all; mean_c_dists];
        clear mean_c_dists
        
        all_block1_dists = [all_block1_dists; block1_dists_mat];
        clear block1_dists_mat
        
    else
    end



    

end
fclose(fid);





%PLOTTING/ANALYSING across datasets
for odor_n = 1:length(odor_list)
    odor_ni = odor_list(odor_n);
    odor_name = dataset(1).stim.odourNames(odor_n).odour;
    odor_resp_mat = saved_sig_areas{1, odor_n};
    curr_plot_color = color_vec(odor_n, :);
        
    %Reality check: does average activity of KC population change across
    %protocol blocks?
    
    fig_handle = figure(1);
    if summary_statistic == 0
        mean_vec = nanmean(odor_resp_mat(:, [3:20] ) );                               %discarding trial 1 for each odor
        se_vec = nanstd(odor_resp_mat(:, [3:20]))./sqrt(size(odor_resp_mat, 1));      %discarding trial 1
        
    elseif summary_statistic == 1
        mean_vec = nanmedian(odor_resp_mat(:, [3:20] ) );                             %discarding trial 1
        se_vec = mad(odor_resp_mat(:, [3:20]))./sqrt(size(odor_resp_mat, 1));      %discarding trial 1
    else
    end
    mean_vec(9:10) = nan;
    se_vec(9:10) = nan;
    %se_vec = nanstd(odor_resp_mat(:, 2:size(odor_resp_mat, 2)))./sqrt(size(odor_resp_mat, 1));      %discarding trial 1
    
    %resizing figure
    
        
    errorbar([1:length(mean_vec)], mean_vec, se_vec, 'Color', curr_plot_color, 'LineWidth', 3)
    hold on
    plot(mean_vec, ['O' curr_plot_color], 'markerfacecolor', curr_plot_color, 'markeredgecolor', curr_plot_color)
    
    set(fig_handle, 'Position', [100, 100, plot_height, plot_width]);
    set(gca,'FontSize', axis_font_size)
    
    xlabel('Trial number')
    if summary_statistic == 0
        ylabel('Mean response area across cells')
    elseif summary_statistic == 1
        ylabel('Median response area across cells')
    else
    end
    
    
    
    %plotting mean and error bar for each cell's responses in first v/s
    %last blocks
    fig_handle2 = figure(2);
    an_tr_list1 = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, 1, blk_trial_window, 0);        %list of trial_window number of trials for this odor from the end of block1 
    an_tr_list2 = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, length(prot_switch_trials), blk_trial_window, 1);        %list of trial_window number of trials for this odor from the end of the last block
    an_tr_list = [an_tr_list1; an_tr_list2];
    
    if summary_statistic == 0
        mean_vec_first = nanmean(odor_resp_mat(:, an_tr_list1), 2);
        mean_vec_last = nanmean(odor_resp_mat(:, an_tr_list2), 2);
        se_vec_first = nanstd(odor_resp_mat(:, an_tr_list1), [], 2)./sqrt(length(an_tr_list));    
        se_vec_last = nanstd(odor_resp_mat(:, an_tr_list2), [], 2)./sqrt(length(an_tr_list));    
    elseif summary_statistic == 1
        mean_vec_first = nanmedian(odor_resp_mat(:, an_tr_list1), 2);
        mean_vec_last = nanmedian(odor_resp_mat(:, an_tr_list2), 2);
        se_vec_first = mad(odor_resp_mat(:, an_tr_list1), [], 2)./sqrt(length(an_tr_list));    
        se_vec_last = mad(odor_resp_mat(:, an_tr_list2), [], 2)./sqrt(length(an_tr_list));    
    else
    end
        
    
    errorbar(mean_vec_first, mean_vec_last, se_vec_first, 'O', 'MarkerFaceColor', curr_plot_color, 'MarkerEdgeColor', curr_plot_color, 'MarkerSize', 3, 'Color', curr_plot_color)
    hold on
    h = herrorbar(mean_vec_first, mean_vec_last, se_vec_last, ['.' curr_plot_color])
    set(h, 'Color', curr_plot_color)
    xlabel('Single cell responses in first block')
    ylabel('Single cell responses in last block')
    %set(gcf, 'Color ', 'w')
    
    %fitting lines to responses in block1 v/s last block
    line_pars = polyfit(mean_vec_first, mean_vec_last, 1);
    linex = [1:1:max([mean_vec_first; mean_vec_last])];
    liney = polyval(line_pars, linex);
    plot(linex, liney, '-', 'Color', curr_plot_color, 'LineWidth', 2)
    
    %plotting slope - 1 reference line for figure 2
    axisi = axis;
    del_max = max(axisi);
    del_min = min(axisi);
    axis([del_min, del_max, del_min, del_max])
    plot([del_min, del_max], [del_min, del_max], '--', 'Color', [0.60, 0.60, 0.60])
    set(fig_handle2, 'Position', [100, 100, plot_height, plot_width]);
    set(gca,'FontSize', axis_font_size)
        

    %STATISTICS (single-cell level analysis)
    %checking if response dists are normal
    %tabulating counts of significant cells
    T = [];

    odor_ni = odor_list(odor_n);
    odor_name = dataset(1).stim.odourNames(odor_ni).odour;
    
    odor_num(odor_n) = odor_n;
    block1_sig_cells(odor_n, 1) = nansum(sig_cell_block_mat(:, odor_ni, 1));
    block3_sig_cells(odor_n, 1) = nansum(sig_cell_block_mat(:, odor_ni, 3));
    del = sig_cell_block_mat(:, :, 3) - sig_cell_block_mat(:, :, 1);
    block1not3_sig_cells(odor_n, 1) = length(find(del(:, odor_ni) == -1) );
    block3not1_sig_cells(odor_n, 1) = length(find(del(:, odor_ni) == 1) );
    block1n3_sig_cells(odor_n, 1) = block1_sig_cells(odor_n) - block1not3_sig_cells(odor_n);
    T_line = [{odor_name}, block1_sig_cells(odor_n, 1), block3_sig_cells(odor_n, 1), block1not3_sig_cells(odor_n, 1), block3not1_sig_cells(odor_n, 1), block1n3_sig_cells(odor_n, 1)];
    T = [T; T_line];
    T = [{'Odor'}, {'Block1'}, {'Block3'}, {'Block1not3'}, {'Block3not1'}, {'Block1and3'}; T];
    disp('summary significant cell counts')
    disp(T)

    data_to_save{1, odor_n} = [mean_vec_first, mean_vec_last];
    
    
    
    
end


%tabulating sig_counts for cells responsive to odors 1 as well as 5 (oCT, MCH)
bothOdSigBlk1 = find(sig_cell_block_mat(:, 1, 1) == 1 & sig_cell_block_mat(:, 5, 1) == 1);
bothOdSigBlk3 = find(sig_cell_block_mat(:, 1, 3) == 1 & sig_cell_block_mat(:, 5, 3) == 1);
bothOdSigBlk1and3 = intersect(bothOdSigBlk1, bothOdSigBlk3);

disp('significant cell counts for OCT and MCH')
T = [length(bothOdSigBlk1), length(bothOdSigBlk3), length(bothOdSigBlk1and3)];
disp([{'Block1'}, {'Block3'}, {'Block1and3'}])
disp(T)


odor1_data = data_to_save{1, 1};
odor2_data = data_to_save{1, 2};
 
% odor1_norm = odor1_data./repmat(median(odor1_data, 1), size(odor1_data, 1), 1);
% odor2_norm = odor2_data./repmat(median(odor2_data, 1), size(odor2_data, 1), 1);

 odor1_norm = odor1_data;       %absolute difference in dF/F, not normalised for each cell
 odor2_norm = odor2_data;       %absolute difference in dF/F, not normalised for each cell


odor1_normdiff = odor1_norm(:, 2) - odor1_norm(:, 1);
odor2_normdiff = odor2_norm(:, 2) - odor2_norm(:, 1);
disp('testing response change distributions')
[h, p] = ttest2(odor1_normdiff, odor2_normdiff)
[h2, p2] = kstest2(odor1_normdiff, odor2_normdiff)

fig_handle3 = figure(3);
[odor_hist(:, 1), hist_centers(:, 1)] = hist(odor1_normdiff, 5);
[odor_hist(:, 2), hist_centers(:, 2)] = hist(odor2_normdiff, 5);
% odor_hist(:, 1) = cumsum(odor_hist(:, 1));
% odor_hist(:, 2) = cumsum(odor_hist(:, 2));

odor_hist(:, 1) = odor_hist(:, 1)./sum(odor_hist(:, 1)).*100;           %normalising numbers of cells in each bin
odor_hist(:, 2) = odor_hist(:, 2)./sum(odor_hist(:, 2)).*100;           %normalising numbers of cells in each bin

for odor_n = 1:n_odors
    plot(hist_centers(:, odor_n), odor_hist(:, odor_n), 'Color', color_vec(odor_n, :), 'LineWidth', 3)       %divided by three bec areas are calculated by summing 25 frames
    hold on
end


set(fig_handle3, 'Position', [100, 100, plot_height, plot_width]);
set(gca,'FontSize', axis_font_size)
xlabel('Change in response size')
ylabel('Percentage of cells')



fig_handle5 = figure(5);
%boxplot(mean_c_dists_all, 'boxstyle', 'filled', 'colors', color_vec)
boxplot(all_c_dists, 'boxstyle', 'filled', 'colors', color_vec, 'widths', 0.2, 'symbol', 'o', 'medianstyle', 'target')
%plotting individual points with random jitter in x
hold on

for odor_n = 1:length(odor_list)
    curr_color = color_vec(odor_n, :);
    rand_vec = (rand(1, size(all_c_dists, 1)).*0.5 + odor_n - 0.25);
    
    plot(rand_vec, all_c_dists, 'o', 'Color', curr_color)
    
end


set(fig_handle5, 'Position', [100, 100, plot_height, plot_width]);
set(gca,'FontSize', axis_font_size)
ylabel('Distances from pre-training centroid')
xlabel('Filler axis label')

%distances of post-pairing trials from pre-pairing centroid
disp('distances of post-pairing trials from pre-pairing centroid')
[h, p] = ttest(all_c_dists(:, 1), all_c_dists(:, 2))

%distances of pre-pairing centroids from origin
%Can't really use this since each distance is measured in a different
%dimensionality. Not normalised because these distances are what are used
%for the normalisation in the first place.
[h, p] = ttest2(centroid_dists_saved(1:(length(centroid_dists_saved)./2), 1), centroid_dists_saved( ((length(centroid_dists_saved)./2)+1):length(centroid_dists_saved), 1));

%distances of post-pairing centroids from origin
centroid_dists_saved(:, 2) = centroid_dists_saved(:, 2)./centroid_dists_saved(:, 1);        %normalising
disp('distances of post-pairing centroids from origin')
[h, p] = ttest2(centroid_dists_saved(1:(length(centroid_dists_saved)./2), 2), centroid_dists_saved( ((length(centroid_dists_saved)./2)+1):length(centroid_dists_saved), 2))



