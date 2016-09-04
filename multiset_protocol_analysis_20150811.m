clear all
close all
list_direc = ['C:\Data\CSHL\dataset_list_Toshi_KC_DN_Led_20150708.txt'];
fid = fopen(list_direc);
direc_counter = 0;

plot_summary_an_only = 0;       %1 - only plots analyses across datasets, 0 - also plots for individual datasets
blk_trial_window = 4;               %uses the last 'trial_window' trials in block1 and the first 'trial_window' trials in blockLast for analysis, plotting. If trial_window is empty, the entire dataset is analysed.
summary_statistic = 1;          % 0 - use mean; 1 - use median

sig_cell_block_mat = [];
color_vec = [[.35, .55, .75]; [.75, .55, .35]; [.55, .35, .75]; [.35, .55, .75]; [.75, .55, .35]; [.75, .35, .55]];

%figure property initialisation variables
plot_height = 300;
plot_width = 300;
axis_font_size = 15;


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
    
%     if isempty(trial_window) == 1
%         trial_window = 1:1:n_trials;
%     else
%         block1_trs = [1, prot_switch_trials(1)];
%         blockend_trs = [(prot_switch_trials(length(prot_switch_trials) - 1) + 1), prot_switch_trials(length(prot_switch_trials))];
%         
%         trial_window = [(block1_trs(2) - trial_window):1:block1_trs(2), blockend_trs(1):1:(blockend_trs(1) + trial_window)];        %list of last 'trial_window' trials of block1 and first 'trial_window' trials of block end
%         
%     end
    
        
    %identifying sig responses on a single trial basis, and then sig
    %responder cells in any individual block
    [resp_areas, sig_trace_mat, sig_cell_mat, sig_cell_block_mati] = cal_sig_responses_20150717(dataset, dff_data_mat, stim_mat, prot_switch_trials);
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
        
        %odor_resps(non_odor_tr_list) = [];
        
        saved_sig_areas{1, odor_n} = [saved_sig_areas{1, odor_n}; odor_resps];
        
        
        
%         if odor_n == 1 && direc_counter == 1
%             sig_cell_traces = zeros(1, 1, 2); 
%         else
%         end
%         
%         
%         sig_cell_traces(:, :, :, odor_n) = saved_sig_traces(:, :, :, odor_n);
%         sig_cell_traces(:, (size(saved_sig_traces, 2) + 1):(size(sig_cell_traces, 2) + length(sig_cell_list)), :) = dff_data_mat(:, sig_cell_list, :, odor_ni);
%         saved_sig_traces(:, :, :, odor_n) = sig_cell_traces;
        
%         keyboard
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

        %Plotting response areas over all trials
        a = colormap('bone');
        a = flipud(a);
        colormap(a)
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

            %plot_data = squeeze(nanmean(dff_all_PCs(:, :, :, odor_n), 3));
            for plot_odor_n = 1:n_odors
                odor_ni = odor_list(plot_odor_n);
                pre_tr_list = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, 1, blk_trial_window, 0);
                post_tr_list = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, length(prot_switch_trials), blk_trial_window, 1);

                plot_data = (dff_all_PCs((stim_frame - 2):(stim_frame + 8), :, pre_tr_list, plot_odor_n));
                trajectory_plotter(plot_data, 4, [], 2, (color_vec(plot_odor_n, :) + 0.15), 0, 0);

                del = input('press enter for post pairing trials');
                plot_data = (dff_all_PCs((stim_frame - 2):(stim_frame + 8), :, post_tr_list, plot_odor_n));
                figure(4)
                trajectory_plotter(plot_data, 4, [], 2, (color_vec(plot_odor_n, :) - 0.2), 0, 0);

            end
            set(gcf, 'Position', [100, 100, plot_height, plot_width]);
            set(gca,'FontSize', axis_font_size)
        
        end

      keyboard
      figure(4)
      close gcf
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
    
    fig_handle = figure(1)
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
    fig_handle2 = figure(2)
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
        
    
    
    errorbar(mean_vec_first, mean_vec_last, se_vec_first, '.', 'Color', curr_plot_color)
    hold on
    h = herrorbar(mean_vec_first, mean_vec_last, se_vec_last, ['.' curr_plot_color])
    set(h, 'Color', curr_plot_color)
    xlabel('Single cell responses in first block')
    ylabel('Single cell responses in last block')
    %set(gcf, 'Color ', 'w')
    
    %fitting lines to responses in block1 v/s last block
    line_pars = polyfit(mean_vec_first, mean_vec_last, 1)
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
    odor_name = dataset(1).stim.odourNames(odor_n).odour;
    
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
    
    
    %calculating a score for distances between response peaks within 
    
    
    
end



odor1_data = data_to_save{1, 1};
odor2_data = data_to_save{1, 2};

% odor1_ratio = odor1_data(:, 1)./odor1_data(:, 2);
% odor2_ratio = odor2_data(:, 1)./odor2_data(:, 2);

odor1_norm = odor1_data./repmat(max(odor1_data, [], 2), 1, 2);
odor2_norm = odor2_data./repmat(max(odor2_data, [], 2), 1, 2);

odor1_normdiff = odor1_norm(:, 2) - odor1_norm(:, 1);
odor2_normdiff = odor2_norm(:, 2) - odor2_norm(:, 1);

[h, p] = ttest2(odor1_normdiff, odor2_normdiff)

fig_handle3 = figure(3)
[odor_hist(:, 1), hist_centers(:, 1)] = hist(odor1_normdiff.*100, 10);
[odor_hist(:, 2), hist_centers(:, 2)] = hist(odor2_normdiff.*100, 10);

for odor_n = 1:n_odors
    plot(hist_centers(:, odor_n), odor_hist(:, odor_n), 'Color', color_vec(odor_n, :), 'LineWidth', 3)
    hold on
    
end


set(fig_handle3, 'Position', [100, 100, plot_height, plot_width]);
set(gca,'FontSize', axis_font_size)
xlabel('Percentage change in response area post-pairing')
ylabel('Number of cells')





