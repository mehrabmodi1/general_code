clear all
close all

dataset_list_paths = [%{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_d5HT1b_Gamma.xls'};...
                      {'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_d5HT1b_AlphaBeta.xls'} ...
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
   
    %loop to go through all experiment datasets listed in list file
    for dir_n = 1:n_dirs
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

        cd(curr_dir);
        tif_name = dir('*.tif');
        stack_obj = ScanImageTiffReader([curr_dir, tif_name(1).name]);
        [frame_time, zoom, n_chans, PMT_offsets] = SI_tif_info(stack_obj);

        %loading extracted raw fluorescence data matrices written by raw_dff_extractor
        raw_data_mat = load([curr_dir 'extracted_raw_data_mat.mat']);
        raw_data_mat = raw_data_mat.raw_data_mat;           %raw F traces extracted from ROIs
        raw_data_mat_orig = raw_data_mat;
        raw_data_mat = raw_data_mat(:, :, stim_mat_simple(:, 1));       %making sure only time-stamp matched trials are used for further analysis
        n_cells = size(raw_data_mat, 2);

        %calculating dF/F traces from raw data
        filt_time = 0.5;            %in s, the time window for boxcar filter for generating filtered traces
        [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, curr_dir);

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
        
        
        %% Analysing pop-representation differences
        for dur_n = 1:length(odor_dur_list)
            curr_dur = odor_dur_list(dur_n);
            summed_sig_cell_vec = sum(sig_cell_mat(:, :, dur_n), 2, 'omitnan');
            all_sig_cells = find(summed_sig_cell_vec > 0);     %all the sig cells across odours, for for the current duration
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
            
            ave_mats = (mean(trace_mats_ds, 3, 'omitnan'));         %normalising mean resps to across-trials sd at each time point
            [del, max_i_mat] = max(ave_mats, [], 1);                       %down-sampled frame numbers of peaks in ave response trace for each cell
            
            %computing across repeat sd of peak responses for each cell
            sd_vec = [];
            for cell_n = 1:length(all_sig_cells)
                for odor_n = 1:n_odors
                    curr_pki = max_i_mat(1, cell_n, odor_n);
                    sd_vec(1, cell_n, odor_n) = std(trace_mats_ds(curr_pki, cell_n, :, odor_n), 1, 3, 'omitnan');
                end
            end
            
            sd_mats = repmat(sd_vec, size(ave_mats, 1), 1, 1);
            ave_mats = ave_mats./sd_mats;
            
            %computing pop-vec distances
            dists_saved = zeros(size(ave_mats, 1), 3);
            for t_step_n = 1:size(ave_mats, 1)
                curr_rep_mat = squeeze(ave_mats(t_step_n, :, :));
                dist_mat = squareform(pdist(curr_rep_mat', 'euclidean'));
                dist_vec = [dist_mat(1, 2), dist_mat(1, 3), dist_mat(2, 3)];
                dists_saved(t_step_n, :) = dist_vec;

            end

            figure(1)
            t_vec = 0:step_size:size(trace_mats, 1);
            t_vec(1) = [];
            plot(t_vec, dists_saved)
            stim_frs = compute_stim_frs(stim_mat, curr_trs(1), frame_time);
            ylabel(['odor-odor pop vec distances, ncells ' num2str(length(all_sig_cells))])
            set_xlabels_time(1, frame_time, 5)
            add_stim_bar(1, stim_frs, [0, 0, 0])
            
            
            
%             figure(2)
%             subplot_tight(1, 2, 1)
%             plot(t_vec, dists_saved)
            
            keyboard
            close figure 1
        end
        
    end
end