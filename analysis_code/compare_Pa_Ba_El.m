clear all
close all

dataset_list_paths = [%{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_d5HT1b_Gamma.xls'};...
                      %{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_c739_AlphaBeta.xls'} ...
                      %{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_d5HT1b_Gamma_low_conc.xls'};...
                      %{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_c739_AlphaBeta_low_conc.xls'};...
                      %{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_c305a_AlphapBetap_low_conc.xls'};...
                      {'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_d5HT1b_Gamma_set2.xls'};...
                      {'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_c739_AlphaBeta_set2.xls'}...
                      ];

suppress_plots = 1;
[del, odor_names] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
a = colormap('bone');
global greymap
greymap = flipud(a);

CoM_vecs = [];
od_resp_tracker_all_saved = [];

for list_n = 1:size(dataset_list_paths, 1)
    curr_dir_list_path = dataset_list_paths{list_n, 1};
    [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
    n_dirs = size(dir_list, 1);
    dataset_list_name = findstr(curr_dir_list_path, 'El_');
    dataset_list_name = curr_dir_list_path((dataset_list_name + 1):(end - 4));
    dataset_list_name(1) = [];
    
    %loop to go through all experiment datasets listed in list file
    n_cells_tot = 0;
    saved_pk_resps_all = [];
    
    saved_int_resps_1s_1 = [];
    saved_int_resps_1s_2 = [];
    saved_int_resps_60s_1 = [];
    saved_int_resps_60s_2 = [];
    
    saved_non_int_resps_1s_1 = [];
    saved_non_int_resps_1s_2 = [];
    saved_non_int_resps_60s_1 = [];
    saved_non_int_resps_60s_2 = [];
    
    sig_resp_mat_all = [];

    od_resp_tracker_all = [];
    
    for dir_n = 1:n_dirs
        saved_an_results.scriptname = mfilename('fullpath');
        curr_dir = [dir_list{dir_n, 1}, '\'];
        tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
        tif_times = tif_times.time_stamps;
        [stim_mat, stim_mat_simple, column_heads, color_vec] = load_params_trains(curr_dir, tif_times);    %reading in trial stimulus parameters after matching time stamps to F traces
        
        %analysing only the first 30 trials
        stim_mat_simple(31:end, :) = [];
        stim_mat(31:end) = [];
        
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

        %reading in number of manually counted cells in FOV
        if exist([curr_dir, 'counted_cells.mat']) == 2
            n_cells_manual = load([curr_dir, 'counted_cells.mat']);
            n_cells_manual = n_cells_manual.ROI_mat;
            n_cells_manual = size(n_cells_manual, 3);
        else
        end
        
        %loading extracted raw fluorescence data matrices written by raw_dff_extractor
        raw_data_mat = load([curr_dir 'extracted_raw_data_mat.mat']);
        raw_data_mat = raw_data_mat.raw_data_mat;           %raw F traces extracted from ROIs
        raw_data_mat_orig = raw_data_mat;
        raw_data_mat = raw_data_mat(:, :, stim_mat_simple(:, 1));       %making sure only time-stamp matched trials are used for further analysis
        n_cells = size(raw_data_mat, 2);

        %calculating dF/F traces from raw data
        filt_time = 1;            %in s, the time window for boxcar filter for generating filtered traces
        [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, curr_dir);
        del = find(dff_data_mat_f < -1);
        dff_data_mat_f(del) = -1;       %forcing crazy values to sane ones
        
        %identifying significantly responsive cells
        [resp_areas, sig_trace_mat, sig_trace_mat_old, sig_cell_mat] = cal_sig_responses_res(dff_data_mat, stim_mat, stim_mat_simple, curr_dir, frame_time);

        %Running data quality control checks
%         sig_cell_mat_old = sig_cell_mat;
%         [sig_cell_mat, all_bad_trs] = cell_data_quality_control(dff_data_mat_f, stim_mat, stim_mat_simple, sig_cell_mat, 1, frame_time);
%         dff_data_mat(:, :, all_bad_trs) = nan;
        %disp([num2str((length(all_bad_trs)./size(dff_data_mat, 3)).*100) ' percent of trials were auto-identified as bad and removed.']);
        sig_cells = find(sum(sum(sig_cell_mat, 3), 2) > 0);         %list of all cells significant for any odor for any duration
        
        %analysing only the first 30 trials
        dff_data_mat(:, :, 31:end) = [];
        
        %saving stuff for sharing data
        if list_n == 1
            share_path = ['C:\Data\Data\Analysed_data\data_sharing\Gamma_set2\fly' int2str(dir_n)];
        elseif list_n == 2
            share_path = ['C:\Data\Data\Analysed_data\data_sharing\AlphaBeta_set2\fly' int2str(dir_n)];
        else
        end
        mkdir(share_path);
        save([share_path, '\dFF_data.mat'], 'dff_data_mat_f');
        save([share_path, '\sig_cells.mat'], 'sig_cell_mat');
        save([share_path, '\stim_mat.mat'], 'stim_mat')
        
        %% Comparing basic stuff
        n_responders = sum(sig_cell_mat(:, odor_list, :), 1, 'omitnan');
        sparsenesses = n_responders./n_cells_manual;
        sparsenesses_saved(:, :, dir_n) = sparsenesses;
                
        od_pair_list = [1, 2; 1, 3; 2, 3];
        for dur_n = 1:length(odor_dur_list)
            curr_dur = odor_dur_list(dur_n);
            
            for odor_pair_n = 1:size(od_pair_list, 1)
                od_pair = od_pair_list(odor_pair_n, :);
                sig_list1 = find(sig_cell_mat(:, odor_list(od_pair(1)), dur_n) == 1);
                sig_list2 = find(sig_cell_mat(:, odor_list(od_pair(2)), dur_n) == 1);
                expected_int = (length(sig_list1)./n_cells).*(length(sig_list2)./n_cells).*n_cells;
                saved_intersections(odor_pair_n, dur_n, dir_n) = length(intersect(sig_list1, sig_list2))./expected_int;
                saved_intersections_n(odor_pair_n, dur_n, dir_n) = length(intersect(sig_list1, sig_list2))./n_cells;
                saved_non_intersections(odor_pair_n, dur_n, dir_n) = (length(union(sig_list1, sig_list2)) - length(intersect(sig_list1, sig_list2)))./n_cells;
                
                
                %saving mean resp traces from intersections/non_intersections
                for pairmem = 1:2                            %only looking at members of the current odour pair
                    odor_np = od_pair_list(odor_pair_n, pairmem);
                    odor_ni = odor_list(odor_np);
                    curr_trs = find(stim_mat_simple(:, 2) == odor_ni & stim_mat_simple(:, 3) == curr_dur);
                    
                    %separately saving responses of each cell to each pairmember
                    if pairmem == 1
                        int_resp_mat1 = mean(dff_data_mat_f(:, intersect(sig_list1, sig_list2), curr_trs), 3, 'omitnan');      
                        non_int_resp_mat1 = mean(dff_data_mat_f(:, setxor(sig_list1, sig_list2), curr_trs), 3, 'omitnan');
                    elseif pairmem == 2
                        int_resp_mat2 = mean(dff_data_mat_f(:, intersect(sig_list1, sig_list2), curr_trs), 3, 'omitnan');
                        non_int_resp_mat2 = mean(dff_data_mat_f(:, setxor(sig_list1, sig_list2), curr_trs), 3, 'omitnan');
                    else
                    end
                end
                
                if dur_n == 1
                    saved_int_resps_1s_1 = [saved_int_resps_1s_1, int_resp_mat1];
                    saved_int_resps_1s_2 = [saved_int_resps_1s_2, int_resp_mat2];
                    saved_non_int_resps_1s_1 = [saved_non_int_resps_1s_1, non_int_resp_mat1];
                    saved_non_int_resps_1s_2 = [saved_non_int_resps_1s_2, non_int_resp_mat2];
                    
                elseif dur_n == 2
                    saved_int_resps_60s_1 = [saved_int_resps_60s_1, int_resp_mat1];
                    saved_int_resps_60s_2 = [saved_int_resps_60s_2, int_resp_mat2];
                    saved_non_int_resps_60s_1 = [saved_non_int_resps_60s_1, non_int_resp_mat1];
                    saved_non_int_resps_60s_2 = [saved_non_int_resps_60s_2, non_int_resp_mat2];
                    

                else
                end
                
            end
            
            
            
            
            saved_pk_resps = [];
            for odor_n = 1:length(odor_list)
                sig_list = find(sig_cell_mat(:, odor_list(odor_n), dur_n) == 1);
                odor_ni = odor_list(odor_n);
                curr_trs = find(stim_mat_simple(:, 2) == odor_ni & stim_mat_simple(:, 3) == curr_dur);
                
                try
                    stim_frs = compute_stim_frs(stim_mat, curr_trs(1), frame_time);
                catch
                    keyboard
                end
                
                ave_mat = mean(dff_data_mat(stim_frs(1):(stim_frs(2) + round(5./frame_time) ), sig_list, curr_trs), 3, 'omitnan');
                pk_resps = max(ave_mat, [], 1, 'omitnan');
                pk_resps = reshape(pk_resps, [], 1);
                saved_pk_resps = [saved_pk_resps; pk_resps];
                clear ave_mat
            end
            
            saved_pk_resps_all = [saved_pk_resps_all; saved_pk_resps];
            
        end
             
        
        %% Doing some Center of Mass based analyses
        for odor_n = 1:length(odor_list)
            odor_ni = odor_list(odor_n);
            dur_n = find(odor_dur_list == 60);
            curr_trs = find(stim_mat_simple(:, 2) == odor_ni & stim_mat_simple(:, 3) == 60);
            curr_sig_cells = find(sig_cell_mat(:, odor_ni, dur_n) == 1);
            
            sig_resp_mat = mean(dff_data_mat_f(:, curr_sig_cells, curr_trs), 3, 'omitnan');
            sig_resp_mat_all = pad_n_concatenate(sig_resp_mat_all, sig_resp_mat, 2, nan);
            
            if odor_ni == 3 || odor_ni == 10
                %checking if current sig cells responded to PA (scored 1), BA (scored 2) or both (scored 3). Other cells are scored 0.
                od_resp_tracker = zeros(size(curr_sig_cells, 1), 1);
                
                PA_respondersi = find(sig_cell_mat(curr_sig_cells, 3, 2) == 1);
                od_resp_tracker(PA_respondersi) = od_resp_tracker(PA_respondersi) + 1;
                
                BA_respondersi = find(sig_cell_mat(curr_sig_cells, 10, 2) == 1);
                od_resp_tracker(BA_respondersi) = od_resp_tracker(BA_respondersi) + 2;
                                
            else
                od_resp_tracker = zeros(size(curr_sig_cells, 1), 1);        %all zeros for EL respondses
            end
            
            od_resp_tracker_all = [od_resp_tracker_all; od_resp_tracker];
        end
        
        
        
        
        %% Analysing pop-representation differences
        %Step0. Computing the maximum, repeat averaged response size for each cell in the long duration stimulus to be used as a
        %normalization factor for each cell's trial averaged, dF/F response in all analyses below.
        for odor_n = 1:length(odor_list)
            odor_ni = odor_list(odor_n);
            curr_trs = find(stim_mat_simple(:, 2) == odor_ni & stim_mat_simple(:, 3) == 60);
            [stim_frs] = compute_stim_frs(stim_mat, curr_trs(1), frame_time); 
            ave_mat(:, :, odor_n) = mean(dff_data_mat(stim_frs(1):(stim_frs(2) + round(10./frame_time)), :, curr_trs), 3, 'omitnan');
            
        end
        max_mat(1, :, :) = max(ave_mat, [], 1);
        max_resp_vec = max(max_mat, [], 3);
        clear ave_mat
        clear max_mat
        
        %Step1. Computing max distance from origin of each n-D odor response vector and then finding the mean across odors to use as a
        %normalization factor for all distances measured in this space.
        dist_type = 'euclidean';            %type of distance to be computed by pdist
        norm_dist_vec = zeros(length(odor_dur_list), 1);
        for dur_n = 2
            curr_dur = odor_dur_list(dur_n);
            summed_sig_cell_vec = sum(sig_cell_mat(:, :, dur_n), 2, 'omitnan');
            all_sig_cells = find(summed_sig_cell_vec > 0);     %all the sig cells across odours, for for the current duration
            if length(all_sig_cells) < 2
                continue
            else
            end
            
            med_dist_vec = zeros(length(odor_list), 1);
            for odor_n = 1:length(odor_list)
                odor_ni = odor_list(odor_n);
                curr_trs = find(stim_mat_simple(:, 2) == odor_ni & stim_mat_simple(:, 3) == curr_dur);
                [stim_frs] = compute_stim_frs(stim_mat, curr_trs(1), frame_time);
                mean_resp_mat = mean(dff_data_mat_f(stim_frs(1):(stim_frs(2) + round(5./frame_time)), all_sig_cells, curr_trs), 3, 'omitnan');
                mean_resp_mat = [zeros(size(mean_resp_mat, 1), 1), mean_resp_mat];
                dist_mat = squareform(pdist(mean_resp_mat', dist_type));
                med_dist_vec(odor_n, 1) = median(dist_mat(2:end, 1));               %median distance of pop-vec from origin across all stim time-points
            end
            norm_dist_vec = max(med_dist_vec);                                      %choosing largest median distance from origin as the normalization distance for this fly
        end
        
        %Step2:Computing various distances and normalising to the pop-vec distance from origin computed in step 1.
        for dur_n = length(odor_dur_list):-1:1
            curr_dur = odor_dur_list(dur_n);
            if dur_n == 2
                summed_sig_cell_vec = sum(sig_cell_mat(:, :, dur_n), 2, 'omitnan');
                all_sig_cells = find(summed_sig_cell_vec > 0);     %all the sig cells across odours, for for the current duration
            else
            end
            %all_sig_cells = 1:1:n_cells;
            if length(all_sig_cells) < 2
                continue
            else
            end
            ave_mats = zeros(size(dff_data_mat, 1), length(all_sig_cells), length(odor_list));
            trace_mats = [];
            for odor_n = 1:length(odor_list)
                odor_ni = odor_list(odor_n);
                curr_trs = find(stim_mat_simple(:, 2) == odor_ni & stim_mat_simple(:, 3) == curr_dur);
                
                %assembling the data matrix for all odors, for current duration
                trace_mats(:, :, :, odor_n) = dff_data_mat_f(:, all_sig_cells, curr_trs);
                
            end
            %averaging step_size frames
           step_size = round(5./frame_time);       %step-size = 5 s
            trace_mats_ds = [];
            for t_step_n = 1:floor(size(trace_mats, 1)./step_size)
                trace_mats_ds(t_step_n, :, :, :) = mean(trace_mats( ((t_step_n - 1).*step_size + 1):(t_step_n.*step_size), :, :, :) , 1, 'omitnan'); 
            end
            
            ave_mats = (mean(trace_mats_ds, 3, 'omitnan'));                %repeat-averaged response matrix
            norm_mat = repmat(max_resp_vec(1, all_sig_cells), size(ave_mats, 1), 1, 1, size(ave_mats, 4));
            ave_mats = ave_mats./norm_mat;                                 %nomalising ave response traces to maximum measured response of each cell
            
            [del, max_i_mat] = max(ave_mats, [], 1);                       %down-sampled frame numbers of peaks in ave response trace for each cell
            
            %computing across repeats, sd of peak responses for each
            %cell
            sd_vec = [];
            for cell_n = 1:length(all_sig_cells)
                for odor_n = 1:n_odors
                    curr_pki = max_i_mat(1, cell_n, odor_n);
                    sd_vec(1, cell_n, odor_n) = std(trace_mats_ds(curr_pki, cell_n, :, odor_n), 1, 3, 'omitnan');
                end
            end
            sd_mats = repmat(sd_vec, size(ave_mats, 1), 1, 1);
            
            %computing pop-vec distances
            dists_saved = zeros(size(ave_mats, 1), 3);
            for t_step_n = 1:size(ave_mats, 1)
                curr_rep_mat = squeeze(ave_mats(t_step_n, :, :));
                dist_mat = squareform(pdist(curr_rep_mat', dist_type));
                dist_mat = dist_mat./norm_dist_vec;                   %normalising all distances measured by factor computed in Step1.
                dist_vec = [dist_mat(1, 2), dist_mat(1, 3), dist_mat(2, 3)];
                dists_saved(t_step_n, :) = dist_vec;

            end

            figure(1)
            t_vec = 0:step_size:size(trace_mats, 1);
            t_vec(1) = [];
            plot(t_vec, dists_saved, 'lineWidth', 2)
            ax_vals = axis;
            ax_vals(4) = 1.2;
            axis(ax_vals);
            if dur_n == 1
                legend('PA-BA', 'PA-EL', 'BA-EL', 'Location', 'northwest')
            elseif dur_n == 2
                legend('PA-BA', 'PA-EL', 'BA-EL', 'Location', 'northeast')
            else
            end

            stim_frs = compute_stim_frs(stim_mat, curr_trs(1), frame_time);
            ylabel(['odor-odor distances, ' num2str(length(all_sig_cells))])
            set_xlabels_time(1, frame_time, 10)
            fig_wrapup(1)
            add_stim_bar(1, stim_frs, [0.5, 0.5, 0.5])
            disp(dataset_list_name)
            
            if suppress_plots == 0
                keyboard
            else
            end
            
            close figure 1
            
            
        end
        n_cells_tot = n_cells_tot + size(dff_data_mat, 2);
    end
    
    
    disp(['n cells ' int2str(n_cells_tot)]);
    saved_an_results.sparsenesses = sparsenesses_saved;
    saved_an_results.sig_intersections = saved_intersections;
    saved_an_results.sig_intersections_n = saved_intersections_n;
    saved_an_results.sig_nonintersections = saved_non_intersections;
    saved_an_results.pk_responses = saved_pk_resps_all;
    saved_an_results.int_resps_1s = cat(3, saved_int_resps_1s_1, saved_int_resps_1s_2);
    saved_an_results.non_int_resps_1s = cat(3, saved_non_int_resps_1s_1, saved_non_int_resps_1s_2);
    saved_an_results.int_resps_60s = cat(3, saved_int_resps_60s_1, saved_int_resps_60s_2);
    saved_an_results.non_int_resps_60s = cat(3, saved_non_int_resps_60s_1, saved_non_int_resps_60s_2);
    
    tr_1s = find(stim_mat_simple(:, 3) == 1);
    tr_60s = find(stim_mat_simple(:, 3) == 60);
    stim_frs_saved(1, :) = compute_stim_frs(stim_mat, tr_1s(1), frame_time);
    stim_frs_saved(2, :) = compute_stim_frs(stim_mat, tr_60s(1), frame_time);
    saved_an_results.stim_frs_saved = stim_frs_saved;
    
    if isempty(findstr(curr_dir_list_path, 'low_conc')) == 0
        save(['C:\Data\Data\Analysed_data\Analysis_results\Yoshi_PaBaEl_lowflux\', dataset_list_name, '_an_results.mat'], 'saved_an_results');
    elseif isempty(findstr(curr_dir_list_path, 'low_conc')) == 1
        save(['C:\Data\Data\Analysed_data\Analysis_results\Yoshi_PaBaEl\', dataset_list_name, '_an_results.mat'], 'saved_an_results');
    else
    end
        
        
    clear saved_an_results
    clear sparsenesses_saved;
    clear saved_intersections;
    clear saved_intersections_n;
    clear saved_non_intersections;
    
    
    %% Doing some CoM based analysis
    CoM_vec = [];
    curr_tr = find(stim_mat_simple(:, 3) == 60);
    stim_frs = compute_stim_frs(stim_mat, curr_tr(1), frame_time);
    for cell_n = 1:size(sig_resp_mat_all, 2)
        curr_vec = sig_resp_mat_all(:, cell_n);
        if sum(isnan(curr_vec)) == length(curr_vec)
            CoM_vec(1, cell_n) = 1;
            continue
        else
        end

        if sign(min(curr_vec)) == -1
            curr_vec = curr_vec + (min(curr_vec).* -1);
        else
        end
        
        CoMs = centerOfMass(curr_vec');
        
        CoM_vec(1, cell_n) = CoMs(1, 2);
    end
    
    sig_resp_mat_sorted = sort_by_vec(sig_resp_mat_all, CoM_vec, 1);
    norm_vec = max(sig_resp_mat_sorted(stim_frs(1, 1):(stim_frs(1, 2) + 100), :), [], 1);
    sig_resp_mat_sorted_n = sig_resp_mat_sorted./repmat(norm_vec, size(sig_resp_mat_sorted, 1), 1);
    
    figure(1)
    imagesc(sig_resp_mat_sorted_n', [0, 1])
    colormap(greymap)
    set_xlabels_time(1, frame_time, 2)
    fig_wrapup(1)
    add_stim_bar(1, stim_frs, [0.65, 0.65, 0.65]);    
    keyboard
    close figure 1
    
    CoM_vecs = pad_n_concatenate(CoM_vecs, CoM_vec, 1, nan);
    od_resp_tracker_all_saved = pad_n_concatenate(od_resp_tracker_all_saved, od_resp_tracker_all, 2, nan);     %keeping track of PA, BA and PA-BA responders for both AB and G cells.
end

%Plotting CoM histograms
[counts_ab, centers_ab] = hist(CoM_vecs(2, :));
counts_g = hist(CoM_vecs(1, :), centers_ab);

%normalising counts
counts_ab = counts_ab./sum(counts_ab);
counts_g = counts_g./sum(counts_g);

figure(1)
plot(centers_ab, counts_ab, 'Color', [0.4, 0.8, 0.3], 'LineWidth', 3)
hold on
plot(centers_ab, counts_g, 'Color', [0.4, 0.5, 0.8], 'LineWidth', 3)
set_xlabels_time_offset(1, frame_time, 10, -25);
fig_wrapup(1)

%plotting histograms of CoMs for unique or both responders for PA and BA
%Gamma KCs
unique_CoMsi = find(od_resp_tracker_all_saved(:, 1) == 1 | od_resp_tracker_all_saved(:, 1) == 2);
unique_CoMs = CoM_vecs(1, unique_CoMsi);
ovlap_CoMsi = find(od_resp_tracker_all_saved(:, 1) == 3);
ovlap_CoMs = CoM_vecs(1, ovlap_CoMsi);

[counts_unique, centers_unique] = hist(unique_CoMs);
counts_ovlap = hist(ovlap_CoMs, centers_unique);

%normalising counts
counts_unique = counts_unique./sum(counts_unique);
counts_ovlap = counts_ovlap./sum(counts_ovlap);
[hG, pG] = ttest2(unique_CoMs, ovlap_CoMs)

figure(2)
plot(centers_unique, counts_unique, 'Color', [0.4, 0.5, 0.8], 'LineWidth', 3)
hold on
plot(centers_unique, counts_ovlap, 'Color', [0.4, 0.5, 0.8].*0.55, 'LineWidth', 3)
title('Gamma KCs')
set_xlabels_time_offset(2, frame_time, 10, -25)
xlabel('CoM time (s)')
ylabel('frac. cell-odour pairs')
fig_wrapup(2)
disp(['n G ovlap = ', int2str(length(ovlap_CoMs))]);
disp(['n G unique = ', int2str(length(unique_CoMs))]);

%AB KCs
unique_CoMsi = find(od_resp_tracker_all_saved(:, 2) == 1 | od_resp_tracker_all_saved(:, 2) == 2);
unique_CoMs = CoM_vecs(2, unique_CoMsi);
ovlap_CoMsi = find(od_resp_tracker_all_saved(:, 2) == 3);
ovlap_CoMs = CoM_vecs(2, ovlap_CoMsi);

[counts_unique, centers_unique] = hist(unique_CoMs);
counts_ovlap = hist(ovlap_CoMs, centers_unique);

%normalising counts
counts_unique = counts_unique./sum(counts_unique);
counts_ovlap = counts_ovlap./sum(counts_ovlap);
[hAB, pAB] = ttest2(unique_CoMs, ovlap_CoMs)

figure(3)
plot(centers_unique, counts_unique, 'Color', [0.4, 0.8, 0.3], 'LineWidth', 3)
hold on
plot(centers_unique, counts_ovlap, 'Color', [0.4, 0.8, 0.3].*0.55, 'LineWidth', 3)
set_xlabels_time_offset(3, frame_time, 10, -25)
title('Alpha/Beta KCs')
xlabel('CoM time (s)')
ylabel('frac. cell-odour pairs')
fig_wrapup(3)
disp(['n A/B ovlap = ', int2str(length(ovlap_CoMs))]);
disp(['n A/B unique = ', int2str(length(unique_CoMs))]);
