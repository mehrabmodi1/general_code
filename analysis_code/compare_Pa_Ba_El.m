clear all
close all

dataset_list_paths = [{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_d5HT1b_Gamma.xls'};...
                      {'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_c739_AlphaBeta.xls'} ...
                      ];

suppress_plots = 1;
[del, odor_names] = xlsread('C:\Data\Code\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
a = colormap('bone');
global greymap
greymap = flipud(a);


for list_n = 1:size(dataset_list_paths, 1)
    curr_dir_list_path = dataset_list_paths{list_n, 1};
    [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
    n_dirs = size(dir_list, 1);
    dataset_list_name = findstr(curr_dir_list_path, '_');
    dataset_list_name = curr_dir_list_path(dataset_list_name(end):(end - 4));
    dataset_list_name(1) = [];
    
    %loop to go through all experiment datasets listed in list file
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
        filt_time = 0.5;            %in s, the time window for boxcar filter for generating filtered traces
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
        
        
        %% Comparing basic stuff
        n_responders = sum(sig_cell_mat(:, odor_list, :), 1, 'omitnan');
        sparsenesses = n_responders./n_cells_manual;
        sparsenesses_saved(:, :, dir_n) = sparsenesses;
                
        od_pair_list = [1, 2; 1, 3; 2, 3];
        for dur_n = 1:length(odor_dur_list)
            
            for odor_pair_n = 1:size(od_pair_list, 1)
                od_pair = od_pair_list(odor_pair_n, :);
                sig_list1 = find(sig_cell_mat(:, odor_list(od_pair(1)), dur_n) == 1);
                sig_list2 = find(sig_cell_mat(:, odor_list(od_pair(2)), dur_n) == 1);
                expected_int = (length(sig_list1)./n_cells).*(length(sig_list2)./n_cells).*n_cells;
                saved_intersections(odor_pair_n, dur_n, dir_n) = length(intersect(sig_list1, sig_list2))./expected_int;
                saved_intersections_n(odor_pair_n, dur_n, dir_n) = length(intersect(sig_list1, sig_list2))./n_cells;
                saved_non_intersections(odor_pair_n, dur_n, dir_n) = (length(union(sig_list1, sig_list2)) - length(intersect(sig_list1, sig_list2)))./n_cells;
               
            end
            
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
        dist_type = 'cityblock';            %type of distance to be computed by pdist
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
        for dur_n = 1:length(odor_dur_list)
            curr_dur = odor_dur_list(dur_n);
            summed_sig_cell_vec = sum(sig_cell_mat(:, :, dur_n), 2, 'omitnan');
            all_sig_cells = find(summed_sig_cell_vec > 0);     %all the sig cells across odours, for for the current duration
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
            ax_vals(4) = 0.15;
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
            
            keyboard
            close figure 1
            
            
        end
        
    end
    
    saved_an_results.sparsenesses = sparsenesses_saved;
    saved_an_results.sig_intersections = saved_intersections;
    saved_an_results.sig_intersections_n = saved_intersections_n;
    saved_an_results.sig_nonintersections = saved_non_intersections;
    
    save(['C:\Data\Data\Analysed_data\Analysis_results\Yoshi_PaBaEl\', dataset_list_name, '_an_results.mat'], 'saved_an_results');
    
    clear saved_an_results
    clear sparsenesses_saved;
    clear saved_intersections;
    clear saved_intersections_n;
    clear saved_non_intersections;
end