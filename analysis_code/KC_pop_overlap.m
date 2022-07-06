clear all
close all

save_path_base = 'C:\Data\Data\Analysed_data\Analysis_results\KC_transition_logregr\';


KC_act_threshold = 0.2;


%Reading in KC response data
KC_types = [{'d5HT1b'}, {'c305a'}, {'c739'}];

paired_od = 1;     %if 1, PA is treated as A, if 2, BA is treated as A

for KC_type_n = 1:3
    KC_type = KC_types{KC_type_n};
    n_flies = length(dir([save_path_base, KC_type])) - 2;
    
    %result log variables
    bin_wt_vecs_flies = [];
    bin_MBON_acts_flies = [];
    threshlin_wt_vecs_flies = [];
    threshlin_MBON_acts_flies = [];
    non_ovlaps_mat_all = [];
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
        bin_thresh = 0.2;
        [bin_wt_vec, bin_MBON_activations] = binary_wts(KC_resp_data, bin_thresh);
        bin_wt_vecs_flies = pad_n_concatenate(bin_wt_vecs_flies, bin_wt_vec, 2, nan);
        bin_MBON_acts_flies = pad_n_concatenate(bin_MBON_acts_flies, bin_MBON_activations, 2, nan);
        
        
        %2. Thresholded linear weights: any cell with a dF/F response > bin_thresh is a responder, wts = norm.(1./thresh_activity)
        bin_thresh = 0.2;
        [threshlin_wt_vec, threshlin_MBON_activations] = threshlin_wts(KC_resp_data, bin_thresh);
        threshlin_wt_vecs_flies = pad_n_concatenate(threshlin_wt_vecs_flies, threshlin_wt_vec, 2, nan);
        threshlin_MBON_acts_flies = pad_n_concatenate(threshlin_MBON_acts_flies, threshlin_MBON_activations, 2, nan);
        
        
        %3. Measuring mean+/- SD to response sizes over repeats of a given
        %odor and measuring overlap over a range of significance cutoffs (allows analog differences to segregate representations)
        n_SD_range = [1.5, 3];
        [non_ovlaps_mat] = pop_ovlaps(KC_resp_data, n_SD_range);
        
        non_ovlaps_mat_all = pad_n_concatenate(non_ovlaps_mat_all, non_ovlaps_mat, 5, nan);
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

    %MBON activations by each odor through these weight vectors
    fig_h = figure('Name', 'binary MBON activations');
    title(KC_type)
    xlabels = [{'Air-A'}, {'Air-A'''}, {'Air-B'}, {'A''-A'}, {'A-A'''}];
    fig_h = scattered_dot_plot_ttest(bin_MBON_acts_flies', fig_h, .6, 1, 4, [0.65, 0.65, 0.65], 1, [], [], xlabels, 2, [0, 0, 0], 2, 0.05, 0, 1, 'force_mean', [], 0);
    hold on
    xlabel('odor')
    ylabel('model MBON activation')
    fig_wrapup(fig_h, [], [100, 120], 0.6);
    
    
    
    %2. Thresholded-linear weight vectors
    %'weights' ie. how many cells were above threshold for Air-A
    fig_h = figure('Name', 'threshlin weight vectors across flies');
    imagesc(threshlin_wt_vecs_flies);
    colorbar;
    title(KC_type)
    colormap('gray')
    ylabel('cell number')
    xlabel('fly number')

    %MBON activations by each odor through these weight vectors
    fig_h = figure('Name', 'threshlin MBON activations');
    fig_h = scattered_dot_plot_ttest(threshlin_MBON_acts_flies', fig_h, .6, 1, 4, [0.65, 0.65, 0.65], 1, [], [], xlabels, 2, [0, 0, 0], 2, 0.05, 0, 1, 'force_mean', [], 0);
    hold on
    xlabel('odor')
    ylabel('model MBON activation')
    fig_wrapup(fig_h, [], [100, 120], 0.6);
    
    
    
    %3. Overlap measures
    n_cells_vec = size(non_ovlaps_mat_all, 1) - sum(isnan(squeeze(non_ovlaps_mat_all(:, 1, 1, 1, :))));
    od_names = xlabels;
    for ref_od_n = 1:size(KC_resp_data, 3)
        fig_h = figure('Name', ['overlaps with ' od_names{ref_od_n}]);
        
        for cutoff_SD_n = 1:length(n_SD_range)
            curr_novlaps = squeeze(non_ovlaps_mat_all(:, :, ref_od_n, cutoff_SD_n, :));
            curr_novlaps = squeeze(sum(curr_novlaps, 1, 'omitnan'));         %summing n overlapping cells
            curr_novlaps = curr_novlaps./repmat(n_cells_vec, size(curr_novlaps, 1), 1);     %computing fraction of cells from absolute counts
            curr_color = [0.7, 0.7, 0.7] * (1./cutoff_SD_n) + 0.2;
            curr_m_color = [0.8, 0.5, 0.4] * (1./cutoff_SD_n) + 0.1;
            fig_h = scattered_dot_plot_ttest(curr_novlaps', fig_h, .6, 1, 2, curr_m_color, 1, [], [], xlabels, 2, curr_m_color, 2, 0.05, 0, 1, 'force_mean', [], 0);
            hold on
            plot([1:size(KC_resp_data,3)]', mean(curr_novlaps, 2, 'omitnan'), 'Color', curr_m_color)
            ylabel([od_names{ref_od_n}, ' non-overlapping fraction']);
           
        end
        fig_wrapup(fig_h, [], [100, 120], 0.6);
        
    end
    
    keyboard
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
    wt_vec = KC_resp_data(:, 1);     %looking only at paired odor responses
    wt_vec(KC_resp_data(:, 1) >= bin_thresh) = 0;
    wt_vec(KC_resp_data(:, 1) < bin_thresh) = 1;
    

    %computing hypothetical MBON activations
    MBON_activations = KC_resp_data'*wt_vec;
end

function [threshlin_wt_vec, threshlin_MBON_activations] = threshlin_wts(KC_resp_data, bin_thresh)
    %averaging input activity across repeats
    KC_resp_data = squeeze(mean(KC_resp_data, 2, 'omitnan'));
    KC_resp_data_orig = KC_resp_data;
    KC_resp_data(KC_resp_data < bin_thresh) = nan;     %thresholding KC resps before using to adjust wts
   
    %constructing weight vector
    threshlin_wt_vec = ones(size(KC_resp_data, 1), 1);
    threshlin_wt_vec = threshlin_wt_vec./KC_resp_data(:, 1);
    threshlin_wt_vec = threshlin_wt_vec./max(threshlin_wt_vec, [], 'omitnan');     %normalizing non-zero weights
    threshlin_wt_vec(isnan(threshlin_wt_vec)) = 1;
    
    %computing hypothetical MBON activations
    threshlin_MBON_activations = KC_resp_data_orig'*threshlin_wt_vec;
    
end


function [non_ovlap_mat_refs_SDcts] = pop_ovlaps(KC_resp_data, n_SD_range)
    mean_resp_mat = squeeze(mean(KC_resp_data, 2, 'omitnan'));
    SD_resp_mat = squeeze(std(KC_resp_data, [], 2, 'omitnan'));
    
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