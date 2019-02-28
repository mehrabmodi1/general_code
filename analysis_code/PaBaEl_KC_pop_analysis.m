clear all
close all

dataset_list_paths = [%{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_d5HT1b_Gamma.xls'};...
                      %{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_c739_AlphaBeta.xls'} ...
                      %{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_d5HT1b_Gamma_low_conc.xls'};...
                      %{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_c739_AlphaBeta_low_conc.xls'};...
                      %{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_c305a_AlphapBetap_low_conc.xls'};...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_Yoshi_PaBaEl_d5HT1b_Gamma_set2.xls'};...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_Yoshi_PaBaEl_c739_AlphaBeta_set2.xls'}...
                      ];

suppress_plots = 0;
[del, odor_names] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
a = colormap('bone');
global greymap
greymap = flipud(a);
script_name = mfilename;

for list_n = 1:size(dataset_list_paths, 1)
    curr_dir_list_path = dataset_list_paths{list_n, 1};
    [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
    n_dirs = size(dir_list, 1);
    dataset_list_name = findstr(curr_dir_list_path, 'El_');
    dataset_list_name = curr_dir_list_path((dataset_list_name + 1):(end - 4));
    dataset_list_name(1) = [];
    
    big_resp_areas = [];
    big_sig_cell_mat = [];
    big_stim_mat_simple = [];
    ovlap_counts = [];
    venn_mats = [];
    
    for dir_n = 1:n_dirs
        saved_an_results.scriptname = mfilename('fullpath');
        curr_dir = [dir_list{dir_n, 1}, '\'];
        tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
        tif_times = tif_times.time_stamps;
        [stim_mat, stim_mat_simple, column_heads, color_vec, g_tr_list] = load_params_trains(curr_dir, tif_times);    %reading in trial stimulus parameters after matching time stamps to F traces
                
        %Reading in experimental parameters
        odor_list = unique(stim_mat_simple(:, 2) );
        n_odors = length(odor_list);
        odor_dur_list = unique(stim_mat_simple(:, 3) );
        n_od_durs = length(odor_dur_list);
        n_trains = max(stim_mat_simple(:, 11));
        saved_an_results.odor_list = odor_list;
        saved_an_results.odor_dur_list = odor_dur_list;
        
        cd(curr_dir);
        tif_name = dir('*.tif');
        stack_obj = ScanImageTiffReader([curr_dir, tif_name(1).name]);
        [frame_time, zoom, n_chans, PMT_offsets] = SI_tif_info(stack_obj);

       
        %loading extracted raw fluorescence data matrices written by raw_dff_extractor
        raw_data_mat = load([curr_dir 'extracted_raw_data_mat.mat']);
        raw_data_mat = raw_data_mat.raw_data_mat;           %raw F traces extracted from ROIs
        raw_data_mat_orig = raw_data_mat;
        raw_data_mat = raw_data_mat(:, :, stim_mat_simple(:, 1));       %making sure only time-stamp matched trials are used for further analysis
        bad_tr_list = 1:1:size(raw_data_mat, 3);
        bad_tr_list(g_tr_list) = [];
        raw_data_mat(:, :, bad_tr_list) = nan;
        n_cells = size(raw_data_mat, 2);
        
        
        %calculating dF/F traces from raw data
        filt_time = 1;            %in s, the time window for boxcar filter for generating filtered traces
        [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, curr_dir);
        del = find(dff_data_mat_f < -1);
        dff_data_mat_f(del) = -1;       %forcing crazy values to sane ones
        %identifying significantly responsive cells
        [resp_areas, sig_trace_mat, sig_trace_mat_old, sig_cell_mat] = cal_sig_responses_res(dff_data_mat, stim_mat, stim_mat_simple, curr_dir, frame_time);

        sig_cells = find(sum(sum(sig_cell_mat, 3), 2) > 0);         %list of all cells significant for any odor for any duration
        
        if suppress_plots == 0
            %visualising odor sorted responses
            [resp_areas_sorted, fig_h, c_vals_mat] = data_quality_resp_matrix(resp_areas, stim_mat_simple, sig_cell_mat, 1, 1);
            pop_c_vals(dir_n, :) = c_vals_mat;
        
        
            %plotting repeat-ave KC pop resp traces to each odor
            KCpop_odor_resp_plotter(dff_data_mat_f, stim_mat_simple, stim_mat, sig_cell_mat, frame_time, 1, script_name);

        else
        end
        %pooling resp data across flies
        
        ave_resp_areas = zeros(size(dff_data_mat, 2), length(odor_list)) + nan;
        
        for odor_n = 1:length(odor_list)
            odor_ni = odor_list(odor_n);
            curr_trs = find(stim_mat_simple(:, 2) == odor_ni);
            ave_resp_areas(1:size(dff_data_mat, 2), odor_n) = mean(resp_areas(:, curr_trs), 2, 'omitnan');
        end
        
        big_resp_areas = cat(1, big_resp_areas, ave_resp_areas);
        big_sig_cell_mat = cat(1, big_sig_cell_mat, sig_cell_mat);
        
        %counting ovlap_fracs of sig responders for each odor pair
        n_responders1 = length(find(sig_cell_mat(:, 3) == 1));
        n_responders2 = length(find(sig_cell_mat(:, 10) == 1));
        n_responders3 = length(find(sig_cell_mat(:, 10) == 1));
        
        n_responders1_2 = length(find(sig_cell_mat(:, 3) == 1 & sig_cell_mat(:, 10) == 1));
        n_responders1_3 = length(find(sig_cell_mat(:, 3) == 1 & sig_cell_mat(:, 11) == 1));
        n_responders2_3 = length(find(sig_cell_mat(:, 10) == 1 & sig_cell_mat(:, 11) == 1));
                
        n_responders1_2_3 = length(find(sig_cell_mat(:, 3) == 1 & sig_cell_mat(:, 10) == 1 & sig_cell_mat(:, 11) == 1));
        
        ovlap_counts(dir_n, 1) = n_responders1_2./(n_responders1 + n_responders2 - n_responders1_2);
        ovlap_counts(dir_n, 2) = n_responders1_3./(n_responders1 + n_responders3 - n_responders1_3);
        ovlap_counts(dir_n, 3) = n_responders2_3./(n_responders2 + n_responders3 - n_responders2_3);
        
        
        %storing venn diagram matrices
%         [del, img] = vennX_MM([n_responders1, n_responders1_2, n_responders2, n_responders2_3, n_responders3, n_responders1_3, n_responders1_2_3]./n_cells, 0.005, [], 0, 1, 0);
%         venn_mats = pad_n_concatenate(venn_mats, img, 3, 0);
        %keyboard
        close figure 1
        close figure 2
        close figure 3
    end
    
    big_resp_areas_all(list_n).big_resp_areas = big_resp_areas;
    big_sig_cell_mat_all(list_n).big_sig_cell_mat = big_sig_cell_mat;
    
    %identifying sig responders to each odor and plotting their mean resp areas for each odor
    %loop to generate plots for resps to each odor
    
    %identifying sig responder cells for each odor
    for odor_n = 1:3
        odor_ni = odor_list(odor_n);
        sig_responders(odor_n).sig_cell_list = find(big_sig_cell_mat(:, odor_ni) == 1);        
    end
    
    %plotting resp size distributions for each cell-set for each odor
    color_vecs =[0.3, 0.3, 0.8; 0.4, 0.5, 0.7; 0.8, 0.3, 0.3]; 
    
    for plot_odor_n = 1:3
        plot_odor_ni = odor_list(plot_odor_n);
        
        for dist_odor_n = 1:3
            
            if suppress_plots == 1
                figure(plot_odor_n)
            elseif suppress_plots == 0
                figure(plot_odor_n + 3)
            else
            end
            curr_resp_cells = sig_responders(dist_odor_n).sig_cell_list;
            curr_resp_areas = big_resp_areas(curr_resp_cells, plot_odor_n);
            bin_vec = 0:0.25:1.5;
            hist_vec = hist(curr_resp_areas, bin_vec);
            hist_vec = hist_vec./sum(hist_vec);         %normalising counts
            plot(bin_vec, hist_vec, 'lineWidth', 3, 'Color', color_vecs(dist_odor_n, :));
            hold on
            axis([0, 1.5, 0, 0.7])
        end
        curr_od_name = odor_names{plot_odor_ni};
       
        title(['responses to ', curr_od_name])
        xlabel('dF/F response size');
        ylabel('norm. counts')
        legend('PA responders', 'BA responders', 'EL responders') 
        fig_wrapup(plot_odor_n, script_name)
        
    end
    
    
    %plotting ovlap counts
    scattered_dot_plot(ovlap_counts, 4, 2, 3, 8, [0.6, 0.6, 0.6], 0, [], [0, 0, 0], [{'PA-BA'}, {'PA-EL'}, {'BA-EL'}], 1, [1, 0.2, 0.2])
    ylabel('overlap fraction');
    fig_wrapup(4, script_name)
    
    [h1, p1] = ttest(ovlap_counts(:, 1), ovlap_counts(:, 2))
    [h2, p2] = ttest(ovlap_counts(:, 1), ovlap_counts(:, 3))
    [h3, p3] = ttest(ovlap_counts(:, 2), ovlap_counts(:, 3))
    
       
    %plotting Venn diagrams
%     figure(5)
%     for dir_n = 1:n_dirs
%         curr_mat = venn_mats(:, :, dir_n);
%         bk_pixi = find(curr_mat == 0);
%         curr_mat(bk_pixi) = 14;
%         a = imagesc(curr_mat);
%         alpha(a, 0.5./dir_n)
%         hold on
%     end
%     colormap(gray)
        
    keyboard
    close figure 1
    close figure 2
    close figure 3
    close figure 4
    close figure 5
end

%making resp-size scatter plot
G_resps = big_resp_areas_all(1).big_resp_areas;
G_resps = G_resps./repmat(G_resps(:, 2), 1, 3);           %normalising each cell's resps to max across odours
G_sig_vec = sum(big_sig_cell_mat_all(1).big_sig_cell_mat, 2);
G_sig_list = find(G_sig_vec > 0);

AB_resps = big_resp_areas_all(2).big_resp_areas;
AB_resps = AB_resps./repmat(AB_resps(:, 2), 1, 3);         %normalising each cell's resps to max across odours
AB_sig_vec = sum(big_sig_cell_mat_all(2).big_sig_cell_mat, 2);
AB_sig_list = find(AB_sig_vec > 0);

G_colour = [0.13, 0.44, 0.72];
AB_colour = [0.13, 0.71, 0.52];
figure(7)
plot(AB_resps(AB_sig_list, 1), AB_resps(AB_sig_list, 3), '.', 'markerSize', 6, 'markerEdgeColor', AB_colour);
hold on
plot(G_resps(G_sig_list, 1), G_resps(G_sig_list, 3), '.', 'markerSize', 6, 'markerEdgeColor', G_colour);
pbaspect([1, 1, 1]);
xlabel('norm PA responses')
ylabel('norm EL responses')
axis([0, 4, 0, 4])

%finding 90% resp slope boundaries
resp_ratios = AB_resps(:, 1)./AB_resps(:, 3);
low_quantile = quantile(resp_ratios, .1);
high_quantile = quantile(resp_ratios, .90);
high_quant_vals = [0, 0; 8.*high_quantile, 8]; 
low_quant_vals = [0, 0; 8.*low_quantile, 8];
plot(high_quant_vals(:, 1), high_quant_vals(:, 2), '-', 'lineWidth', 2, 'Color', AB_colour )
plot(low_quant_vals(:, 1), low_quant_vals(:, 2), '-', 'lineWidth', 2, 'Color', AB_colour )

resp_ratios = G_resps(:, 1)./G_resps(:, 3);
low_quantile = quantile(resp_ratios, .1);
high_quantile = quantile(resp_ratios, .90);
high_quant_vals = [0, 0; 8.*high_quantile, 8]; 
low_quant_vals = [0, 0; 8.*low_quantile, 8];
plot(high_quant_vals(:, 1), high_quant_vals(:, 2), '-', 'lineWidth', 2, 'Color', G_colour )
plot(low_quant_vals(:, 1), low_quant_vals(:, 2), '-', 'lineWidth', 2, 'Color', G_colour )
fig_wrapup(7, script_name)


%making resp-size scatter plot 2
G_resps = big_resp_areas_all(1).big_resp_areas;
G_resps = G_resps./repmat(G_resps(:, 3), 1, 3);           %normalising each cell's resps to max across odours
G_sig_vec = sum(big_sig_cell_mat_all(1).big_sig_cell_mat, 2);
G_sig_list = find(G_sig_vec > 0);

AB_resps = big_resp_areas_all(2).big_resp_areas;
AB_resps = AB_resps./repmat(AB_resps(:, 3), 1, 3);         %normalising each cell's resps to max across odours
AB_sig_vec = sum(big_sig_cell_mat_all(2).big_sig_cell_mat, 2);
AB_sig_list = find(AB_sig_vec > 0);


figure(8)
plot(AB_resps(AB_sig_list, 1), AB_resps(AB_sig_list, 2), '.', 'markerSize', 6, 'markerEdgeColor', AB_colour);
hold on
plot(G_resps(G_sig_list, 1), G_resps(G_sig_list, 2), '.', 'markerSize', 6, 'markerEdgeColor', G_colour);
pbaspect([1, 1, 1]);
xlabel('norm PA responses')
ylabel('norm BA responses')
axis([0, 4, 0, 4])

%finding 90% resp slope boundaries
resp_ratios = AB_resps(:, 1)./AB_resps(:, 2);
low_quantile = quantile(resp_ratios, .1);
high_quantile = quantile(resp_ratios, .90);
high_quant_vals = [0, 0; 8.*high_quantile, 8]; 
low_quant_vals = [0, 0; 8.*low_quantile, 8];
plot(high_quant_vals(:, 1), high_quant_vals(:, 2), '-', 'lineWidth', 2, 'Color', AB_colour )
plot(low_quant_vals(:, 1), low_quant_vals(:, 2), '-', 'lineWidth', 2, 'Color', AB_colour )

resp_ratios = G_resps(:, 1)./G_resps(:, 2);
low_quantile = quantile(resp_ratios, .1);
high_quantile = quantile(resp_ratios, .90);
high_quant_vals = [0, 0; 8.*high_quantile, 8]; 
low_quant_vals = [0, 0; 8.*low_quantile, 8];
plot(high_quant_vals(:, 1), high_quant_vals(:, 2), '-', 'lineWidth', 2, 'Color', G_colour )
plot(low_quant_vals(:, 1), low_quant_vals(:, 2), '-', 'lineWidth', 2, 'Color', G_colour )
fig_wrapup(8, script_name)




