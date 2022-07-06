clear all
close all

save_path_base = 'C:\Data\Data\Analysed_data\Analysis_results\KC_transition_logregr\';


KC_act_threshold = 0.2;

fig_save_path = 'C:\Users\Mehrab\Dropbox (HHMI)\data_sharing\figure_for_Glenn_20220314\KC_analysis_figs\linear_models_nonovlap\';
save_figs = 0;

%Reading in KC response data
KC_types = [{'d5HT1b'}, {'c305a'}, {'c739'}];

paired_od = 1;     %if 1, PA is treated as A, if 2, BA is treated as A

for KC_type_n = 1:3
    fig_n = 0;
    KC_type = KC_types{KC_type_n};
    n_flies = length(dir([save_path_base, KC_type])) - 2;
    
    %result log variables
    bin_wt_vecs_flies = [];
    bin_MBON_acts_flies = [];
    threshlin_wt_vecs_flies = [];
    threshlin_MBON_acts_flies = [];
    non_ovlaps_mat_all = [];
    KC_resps_plot_all = [];
    KC_wts_plot_all = [];
    KC_binwts_plot_all = [];
    baseline_MBON_acts_flies = [];
    resp_frac_tables_all = [];
    for fly_n = 1:n_flies
        
        %skipping flies with few cells
        if KC_type_n <=2 && fly_n == 4  
            continue
        else
        end
        
        KC_data = load([save_path_base, KC_type, '\fly', num2str(fly_n), '\fly_data.mat']);
        KC_data = KC_data.fly_data;
        KC_resp_data = KC_data.resp_sizes;
        KC_resp_data(KC_resp_data < -0.5) = -0.5;
        
        %correcting ordering of transition odor stimuli based on KC_data.stim_mat_key 
        KC_resp_data_orig = KC_resp_data;
        KC_resp_data(:, :, 4) = KC_resp_data_orig(:, :, 5);
        KC_resp_data(:, :, 5) = KC_resp_data_orig(:, :, 4);       
        
        %re-arranging to use PA and BA as A separately
        KC_resp_data1 = KC_resp_data;   %using PA as A
        KC_resp_data2 = KC_resp_data;
        KC_resp_data2(:, :, 1) = KC_resp_data(:, :, 2);
        KC_resp_data2(:, :, 2) = KC_resp_data(:, :, 1);
        KC_resp_data2(:, :, 4) = KC_resp_data(:, :, 5);
        KC_resp_data2(:, :, 5) = KC_resp_data(:, :, 4);
        
        %switching to BA being used as A if specified
        if paired_od == 2
            KC_resp_data = KC_resp_data2;
        else
        end
        
        %computing population response overlaps based on different rules
        %1. any cell with a dF/F response > bin_thresh is a responder
        bin_thresh = KC_act_threshold;
        [bin_wt_vec, bin_MBON_activations] = binary_wts(KC_resp_data, bin_thresh);
        bin_wt_vecs_flies = pad_n_concatenate(bin_wt_vecs_flies, bin_wt_vec, 2, nan);
        bin_MBON_acts_flies = pad_n_concatenate(bin_MBON_acts_flies, bin_MBON_activations, 2, nan);
        KC_binwts_plot_all = [KC_binwts_plot_all; bin_wt_vec];    %logging wts across flies
        
               
        %2. Thresholded linear weights: any cell with a dF/F response > bin_thresh is a responder, wts = norm.(1./thresh_activity)
        %bin_thresh = 0.2;
        [threshlin_wt_vec, threshlin_MBON_activations, KC_A_resps, baseline_MBON_acts] = threshlin_wts(KC_resp_data, bin_thresh);
        threshlin_wt_vecs_flies = pad_n_concatenate(threshlin_wt_vecs_flies, threshlin_wt_vec, 2, nan);
        threshlin_MBON_acts_flies = pad_n_concatenate(threshlin_MBON_acts_flies, threshlin_MBON_activations, 2, nan);
        KC_resps_plot_all = [KC_resps_plot_all; KC_A_resps];    %logging responses to Air-A across flies
        KC_wts_plot_all = [KC_wts_plot_all; threshlin_wt_vec];    %logging wts across flies
        baseline_MBON_acts_flies = pad_n_concatenate(baseline_MBON_acts_flies, baseline_MBON_acts, 2, nan);     %logging normalized pop activity across flies

             
        %3. Measuring mean+/- SD to response sizes over repeats of a given
        %odor and measuring overlap over a range of significance cutoffs (allows analog differences to segregate representations)
        n_SD_range = [1.5, 20];
        use_only_sigcells = 1;      %1 - use only arbitrarily defined significant responders 0 - use all cells to count non-overlap fraction
        [non_ovlaps_mat] = pop_ovlaps(KC_resp_data, KC_data.sig_cell_mat, n_SD_range, use_only_sigcells);
        
        non_ovlaps_mat_all = pad_n_concatenate(non_ovlaps_mat_all, non_ovlaps_mat, 5, nan);
        
        %compiling responder fractions for each odor and odor pair
        sig_cell_mat = KC_data.sig_cell_mat;
        for od_n1 = 1:size(sig_cell_mat, 2)
            for od_n2 = 1:size(sig_cell_mat, 2)
                n_bothresponders = length(find(sum([sig_cell_mat(:, od_n1), sig_cell_mat(:, od_n2)], 2, 'omitnan') == 2));
                resp_frac_table(od_n1, od_n2) = n_bothresponders;
            end
        end
        resp_frac_table = resp_frac_table./size(sig_cell_mat, 1);
        
        resp_frac_tables_all = pad_n_concatenate(resp_frac_tables_all, resp_frac_table, 3, nan);
    end
    
    %plotting results
    
    %1. Binary weight vectors 
    %'weights' ie. how many cells were above threshold for Air-A
    fig_h = figure('Name', 'weight vectors across flies');
    imagesc(bin_wt_vecs_flies);
    colorbar;
    colormap('gray')
    ylabel('cell number')
    xlabel('fly number')
    
    fig_n = fig_n + 1;
    fig_name = fig_h.Name;
    
    if save_figs == 1
        savefig([fig_save_path, KC_type, ' ', num2str(fig_n), ' ', fig_name, '.fig']);
    else
    end
    
    %MBON activations by each odor through these weight vectors
    fig_h = figure('Name', 'binary MBON activations');
    title(KC_type)
    xlabels = [{'A'}, {'A'}, {'A'''}, {'A'''}, {'B'}, {'B'}, {'A''-A'}, {'A''-A'}, {'A-A'''}, {'A-A'''}];
    %interleaving baseline MBON acts with weight-adjusted MBON acts
    mean_color = [.8, 0.5, 0.3];
    data_mat = [];
    color_vec = [];
    for od_n = 1:size(baseline_MBON_acts_flies, 1)
        data_mat = [data_mat, baseline_MBON_acts_flies(od_n, :)', bin_MBON_acts_flies(od_n, :)'];
        color_vec = [color_vec; [0.7, 0.7, 0.7]; [0, 0, 0]];
    end
    data_mat = data_mat./max(mean(baseline_MBON_acts_flies, 2, 'omitnan'));     %normalizing to largest mean baseline response
    fig_h = scattered_dot_plot_ttest(data_mat, fig_h, .6, [3, 1.5], 4, color_vec, 1, [], [], xlabels, 2, color_vec, 2, 0.05, 0, 1, 'force_mean', [], 0, 5.6);
%     hold on
%     %plotting KC resp data with all weights = 1 to model pre-pairing baseline
%     fig_h = scattered_dot_plot_ttest(baseline_MBON_acts_flies', fig_h, .6, 1, 4, [0.65, 0.65, 0.65], 1, [], [], xlabels, 2, [0.65, 0.65, 0.65], 2, 0.05, 0, 1, 'force_mean', [], 0);
    xlabel('odor')
    ylabel('model MBON activation')
    fig_wrapup(fig_h, [], [60, 80], 0.6);
    fig_n = fig_n + 1;
    fig_name = fig_h.Name;
    
    %Statistical testing
    [p_simp, h] = signrank(data_mat(:, 2), data_mat(:, 4));
    [p_trans, h] = signrank(data_mat(:, 8), data_mat(:, 10));
    
    pvals_bin = bonf_holm([p_simp, p_trans], 0.01)
    
    if save_figs == 1
        savefig([fig_save_path, KC_type, ' ', num2str(fig_n), ' ', fig_name, '.fig']);
    else
    end
    
    %2. Thresholded-linear weight vectors
    %'weights' ie. how many cells were above threshold for Air-A
    fig_h = figure('Name', 'threshlin weight vectors across flies');
    imagesc(threshlin_wt_vecs_flies);
    colorbar;
    title(KC_type)
    colormap('gray')
    ylabel('cell number')
    xlabel('fly number')
    fig_n = fig_n + 1;
    fig_name = fig_h.Name;
    
    %3. Thresholded-linear weights v/s ativity
    %'weights' ie. how many cells were above threshold for Air-A
    fig_h = figure('Name', 'threshlin weight vectors across flies');
    imagesc(threshlin_wt_vecs_flies);
    colorbar;
    title(KC_type)
    colormap('gray')
    ylabel('cell number')
    xlabel('fly number')
    fig_n = fig_n + 1;
    fig_name = fig_h.Name;
    
    if save_figs == 1
        savefig([fig_save_path, KC_type, ' ', num2str(fig_n), ' ', fig_name, '.fig']);
    else
    end
    
    %MBON activations by each odor through these weight vectors
    fig_h = figure('Name', 'threshlin MBON activations');
    data_mat = [];
    color_vec = [];
    for od_n = 1:size(baseline_MBON_acts_flies, 1)
        data_mat = [data_mat, baseline_MBON_acts_flies(od_n, :)', threshlin_MBON_acts_flies(od_n, :)'];
        color_vec = [color_vec; [0.7, 0.7, 0.7]; [0, 0, 0]];
    end
    data_mat = data_mat./max(mean(baseline_MBON_acts_flies, 2, 'omitnan'));     %normalizing to largest mean baseline response
    fig_h = scattered_dot_plot_ttest(data_mat, fig_h, .6, [3, 1.5], 4, color_vec, 1, [], [], xlabels, 2, color_vec, 2, 0.05, 0, 1, 'force_mean', [], 0, 5.6);
%     hold on
%     %plotting KC resp data with all weights = 1 to model pre-pairing baseline
%     fig_h = scattered_dot_plot_ttest(baseline_MBON_acts_flies', fig_h, .6, 1, 4, [0.65, 0.65, 0.65], 1, [], [], xlabels, 2, [0.65, 0.65, 0.65], 2, 0.05, 0, 1, 'force_mean', [], 0);
    xlabel('odor');
    ylabel('model MBON activation');
    fig_wrapup(fig_h, [], [60, 80], 0.6);
    fig_n = fig_n + 1;
    fig_name = fig_h.Name;
    
    %Statistical testing
    [p_simp, h] = signrank(data_mat(:, 2), data_mat(:, 4));
    [p_trans, h] = signrank(data_mat(:, 8), data_mat(:, 10));
    
    pvals_lin = bonf_holm([p_simp, p_trans], 0.01)
    
    if save_figs == 1
        savefig([fig_save_path, KC_type, ' ', num2str(fig_n), ' ', fig_name, '.fig']);
    else 
    end
    
    
    %3. Overlap measures
    n_cells_vec = size(non_ovlaps_mat_all, 1) - sum(isnan(squeeze(non_ovlaps_mat_all(:, 1, 1, 1, :))));
    xlabels = [{'A'}, {'A'''}, {'B'}, {'A''-A'}, {'A-A'''}];
    od_names = xlabels;
    for ref_od_n = 1:size(KC_resp_data, 3)
        fig_h = figure('Name', ['non-overlaps with ' od_names{ref_od_n}]);
        
        for cutoff_SD_n = 1:length(n_SD_range)
            curr_novlaps = squeeze(non_ovlaps_mat_all(:, :, ref_od_n, cutoff_SD_n, :));
            curr_novlaps = squeeze(sum(curr_novlaps, 1, 'omitnan'));         %summing n overlapping cells
            curr_novlaps = curr_novlaps./repmat(n_cells_vec, size(curr_novlaps, 1), 1);     %computing fraction of cells from absolute counts
            curr_color = [0.7, 0.7, 0.7] * (1./cutoff_SD_n) + 0.2;
            curr_m_color = [0.7, 0.7, 0.7] * (1./cutoff_SD_n) + 0.1;
            [fig_h, del, saved_col_centers] = scattered_dot_plot_ttest(curr_novlaps', fig_h, .6, 1, 4, curr_m_color, 1, [], [], xlabels, 2, curr_m_color, 2, 0.05, 0, 1, 'force_mean', [], 0, 5.6);
            hold on
            plot(saved_col_centers, mean(curr_novlaps, 2, 'omitnan'), 'Color', curr_m_color)
            ylabel(['non-overlap with ', od_names{ref_od_n}]);
           
        end
        legend({[num2str(n_SD_range(1)), ' SDs'], [num2str(n_SD_range(2)), ' SDs']});
        fig_wrapup(fig_h, [], [60, 120], 0.6);
        fig_n = fig_n + 1;
        fig_name = fig_h.Name;
        
        if save_figs == 1
            savefig([fig_save_path, KC_type, ' ', num2str(fig_n), ' ', fig_name, '.fig']);
        else
        end
    end
    
    KC_resps_plot_all(KC_resps_plot_all < 0) = 0;
    %4. Plotting binary wts v/s responses to Air-A
    fig_h = figure('Name', 'A-resps v/s bin wts');
    plot(KC_resps_plot_all(:, 1), ones(length(KC_binwts_plot_all), 1), '.',  'Color', [0.7, 0.7, 0.7], 'MarkerSize', 6);
    hold on
    plot(squeeze((KC_resps_plot_all(:, 1))), KC_binwts_plot_all, '.', 'Color', [0, 0, 0], 'MarkerSize', 12);
    xlabel('norm. KC response (dF/F)');
    ylabel('KC weight');
    fig_wrapup(fig_h, [], [5, 7], 0.6);
    
    %5. Plotting threshlin wts v/s responses to Air-A
    fig_h = figure('Name', 'A-resps v/s wts');
    plot(KC_resps_plot_all(:, 1), ones(length(KC_binwts_plot_all), 1), '.',  'Color', [0.7, 0.7, 0.7], 'MarkerSize', 6);
    hold on
    plot(squeeze((KC_resps_plot_all(:, 1))), KC_wts_plot_all, '.', 'Color', [0, 0, 0], 'MarkerSize', 12);
    xlabel('norm. KC response (dF/F)');
    ylabel('KC weight');
    fig_wrapup(fig_h, [], [5, 7], 0.6);
    ax_vals = axis;
    ax_vals(3) = 0;
    axis(ax_vals);
    
    
    
    %compiling responder fractions for each odor and odor pair and
    %outputting a table
    mean_table = mean(resp_frac_tables_all, 3, 'omitnan');
    sd_table = std(resp_frac_tables_all, [], 3, 'omitnan');
    
    %writing tables to file
    table_path = 'C:\Data\Data\Analysed_data\data_sharing\PABAEL_KC_responder_fracs\';
    writematrix(mean_table, [table_path, KC_type, '_resp_frac_means.xlsx']);
    writematrix(sd_table, [table_path, KC_type, '_resp_frac_SDs.xlsx']);
    
    keyboard
    clear resp_frac_tables_all
    close all
end




%worker functions
%1. Binarised weights
function [wt_vec, MBON_activations] = binary_wts(KC_resp_data, bin_thresh)
%this function generates a binary weight vector with zeros for paired odor responders
%larger than bin_thresh

    %averaging input activity across repeats
    KC_resp_data = squeeze(mean(KC_resp_data, 2, 'omitnan'));
   
    %constructing a binary weight vector
    KC_resp_data = KC_resp_data./max(squeeze(KC_resp_data(:, 1)));  %normalizing before threholding
    wt_vec = KC_resp_data(:, 1);     %looking only at paired odor responses
    wt_vec(KC_resp_data(:, 1) >= bin_thresh) = 0;
    wt_vec(KC_resp_data(:, 1) < bin_thresh) = 1;
    

    %computing hypothetical MBON activations
    MBON_activations = (KC_resp_data'*wt_vec)./length(wt_vec);  %dividing by n_cells after a matrix multiplication
end

function [threshlin_wt_vec, threshlin_MBON_activations, KC_A_resps, baseline_MBON_acts] = threshlin_wts(KC_resp_data, bin_thresh)
    %averaging input activity across repeats
    KC_resp_data = squeeze(mean(KC_resp_data, 2, 'omitnan'));
    KC_resp_data_orig = KC_resp_data;
    KC_resp_data = KC_resp_data./max(squeeze(KC_resp_data(:, 1)));  %normalizing before threholding
    KC_A_resps = KC_resp_data(:, 1);        %logging responses for plots
    KC_all_resps = KC_resp_data;
    KC_resp_data(KC_resp_data < bin_thresh) = nan;     %thresholding KC resps before using to adjust wts
    %KC_resp_data = KC_resp_data./max(squeeze(KC_resp_data(:, 1))); %normalizing after threholding
    
    %constructing weight vector
    threshlin_wt_vec = ones(size(KC_resp_data, 1), 1);
    threshlin_wt_vec = threshlin_wt_vec - KC_resp_data(:, 1);   %activity subtracted from weights
    %threshlin_wt_vec = threshlin_wt_vec./KC_resp_data(:, 1);   %weights inversely proportional to activity ie. exponentially reducing with activity
    
    threshlin_wt_vec = threshlin_wt_vec./max(threshlin_wt_vec, [], 'omitnan');     %normalizing weights
    threshlin_wt_vec(isnan(threshlin_wt_vec)) = 1;
    
    %computing hypothetical MBON activations
    threshlin_MBON_activations = (KC_resp_data_orig'*threshlin_wt_vec)./length(threshlin_wt_vec);    %matrix multiplication, not pair-wise
    ones_wt_vec = ones(size(threshlin_wt_vec, 1), 1);
    baseline_MBON_acts = (KC_resp_data_orig'*ones_wt_vec)./length(ones_wt_vec);    %dividing by n_cells after matrix mult
    
    %KC_resp_data(isnan(KC_resp_data)) = 0.000001;
    %KC_A_resps = KC_resp_data(:, 1);
    
    %temporarily plotting wts and activities across flies for later plotting
%     fig_h = figure();
%     KC_resp_data(isnan(KC_resp_data)) = 0.000001;
%     plot(squeeze((KC_resp_data(:, 1))), threshlin_wt_vec, '.', 'Color', [0.7, 0.7, 0.7], 'MarkerSize', 12);
%     xlabel('norm. KC response (dF/F)');
%     ylabel('KC weight');
%     fig_wrapup(fig_h, [], [60, 80], 0.6);
%     keyboard
end


function [non_ovlap_mat_refs_SDcts] = pop_ovlaps(KC_resp_data, sig_cell_mat, n_SD_range, use_only_sigcells)
    


    mean_resp_mat = squeeze(mean(KC_resp_data, 2, 'omitnan'));
    SD_resp_mat = squeeze(std(KC_resp_data, [], 2, 'omitnan'));
    
    if use_only_sigcells == 1
        mean_resp_mat = mean_resp_mat.*sig_cell_mat;
        SD_resp_mat = SD_resp_mat.*sig_cell_mat;
    else
    end
    
    %repeating for each n SD cutoff specified
    non_ovlap_mat_refs_SDcts = [];
    %dim1 - cells, dim2 - odors, dim3 - ref_odors, dim4 - n_SDs
    for SD_count_n = 1:length(n_SD_range)
        n_SDs = n_SD_range(SD_count_n);
        
        %looking for overlap with each odor as reference
        n_ods = size(KC_resp_data, 3);
        non_ovlap_mat_allrefs = [];
        for od_n = 1:n_ods
            curr_mean_dists = abs(mean_resp_mat - repmat(mean_resp_mat(:, od_n), 1, n_ods));    %distance of each cell's mean response to each odor from its mean response to the reference odor
            SD_thresh_mat = SD_resp_mat.*n_SDs; %+ repmat(SD_resp_mat(:, od_n).*n_SDs, 1, n_ods);  %computing summed, SD cutoffs for each cell, for each odor wrt the reference odor     
            
            %defining overlap with reference odor for each cell if mean distance is less than specified n SDs
            ovlap_mat = curr_mean_dists - SD_thresh_mat;        %subtracting cutoff distance from absolute distance between means
            nonovlap_mat = sign(ovlap_mat);
            nonovlap_mat(nonovlap_mat < 0) = 0;           %binary matrix of cell responses not overlapping with reference odor response
            non_ovlap_mat_allrefs = pad_n_concatenate(non_ovlap_mat_allrefs, nonovlap_mat, 3, nan);
        end
        non_ovlap_mat_refs_SDcts = pad_n_concatenate(non_ovlap_mat_refs_SDcts, non_ovlap_mat_allrefs, 4, nan);
    end

end
