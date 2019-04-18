clear all
close all

direc_lists_mat =  [...
                      {'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_fluc_stim_somas_20171226.xls'}...
                        %{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_fluc_stim_axons_20180117'}, ...
                   ]; 

save_path = 'C:\Data\Data\Analysed_data\Analysis_results\train_clustering\';
n_direc_lists = size(direc_lists_mat, 1);
                
global color_vec;                
color_vec = load('C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Code\std_color_vec.txt');
a = colormap('bone');
global greymap
greymap = flipud(a);
colormap(greymap)
suppress_plots = 0;       %0 - doesn't plot quality control stuff, 1 - plots stuff

kernel_width = 5;         %in s, the duration of dF/F trace to use and extract as the Ca-response kernel.
n_trs_to_analyse = 15;

[del, odor_names] = xlsread('C:\Data\Code\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);




%loop to go through all directory lists
for direc_list_n = 1:n_direc_lists

    list_direc = direc_lists_mat{direc_list_n, 1};
    [del, curr_direc_list] = xlsread(list_direc, 1);
    n_dirs = size(curr_direc_list, 1);
    direc_counter = 0;
    
    %parsing direc list path for name of direc list
    namei = findstr(list_direc, '\');
    namei = namei(end) + 1;
    dir_list_name = (list_direc(namei:(end-4)));
    
    MSE_mat = [];
    %loop to go through all experiment datasets listed in list file
    for direc_counter = 1:n_dirs
        %% House-keeping
        direc = curr_direc_list{direc_counter, 1};
        
        direc = [direc, '\'];
        
        %loading in and parsing params file to get stimulus parameter
        %details
        tif_times = load([direc, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
        tif_times = tif_times.time_stamps;
        [stim_mat, stim_mat_simple, column_heads] = load_params_trains(direc, tif_times);
       
        odor_list = unique(stim_mat_simple(:, 2) );
        n_odors = length(odor_list);
        odor_dur_list = unique(stim_mat_simple(:, 3) );
        n_od_durs = length(odor_dur_list);
        n_trains = max(stim_mat_simple(:, 11));
        
        
        %reading ScanImage meta data from a raw .tif
        old_direc = pwd;
        cd(direc);
        tif_names = dir('*.tif');
        cd(old_direc);
        clear old_direc
        tif_name = tif_names(1).name;
        stack_obj = ScanImageTiffReader([direc, tif_name]);
        %metadata = stack_obj.descriptions;
        metadata = stack_obj.metadata;
        fr_ratei = strfind(metadata, 'scanFrameRate');
        frame_rate = metadata((fr_ratei + 16):(fr_ratei + 20));
        frame_rate = str2num(frame_rate);                       %base frame rate in Hz
        avg_factori = strfind(metadata, 'logAverageFactor');
        avg_factor = metadata((avg_factori + 19):(avg_factori + 20));
        avg_factor = str2num(avg_factor);
        frame_rate = frame_rate./avg_factor;
        global frame_time;
        frame_time = 1./frame_rate.*1000;     %in ms
        stim_time = stim_mat_simple(1, 7);
        stim_fr = round((stim_time.*1000)./frame_time);
        post_od_scan_dur = stim_mat_simple(1, 10);
        
        %loading extracted raw fluorescence data matrices written by
        %raw_dff_extractor
        raw_data_mat = load([direc 'extracted_raw_data_mat.mat']);
        raw_data_mat = raw_data_mat.raw_data_mat;           %raw F traces extracted from ROIs
        raw_data_mat_orig = raw_data_mat;
        raw_data_mat = raw_data_mat(:, :, stim_mat_simple(:, 1));       %making sure only time-stamp matched trials are used for further analysis
        
        %calculating dF/F traces from raw data
        filt_time = 200;            %in ms, the time window for boxcar filter for generating filtered traces
        [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, direc);
        
        %identifying significantly responsive cells
        [resp_areas, sig_trace_mat, sig_trace_mat_old, sig_cell_mat] = cal_sig_responses_res(dff_data_mat, stim_mat, stim_mat_simple, direc, frame_time);
        
        %Running data quality control checks
        sig_cell_mat_old = sig_cell_mat;
        [sig_cell_mat, all_bad_trs] = cell_data_quality_control(dff_data_mat_f, stim_mat, stim_mat_simple, sig_cell_mat, suppress_plots);
        disp([num2str((length(all_bad_trs)./size(dff_data_mat, 3)).*100) ' percent of trials were auto-identified as bad and removed.']);
        dff_data_mat(:, :, all_bad_trs) = nan;
        
        %% Running clustering algorithm
        for odor_n = 1:n_odors
            odor_ni = odor_list(odor_n);
            curr_sig_cells = find(sig_cell_mat(:, odor_ni) == 1);
            for train_n = 1:n_trains
                other_train_n = 1:n_trains;
                del = find(other_train_n == train_n);
                other_train_n(del) = [];
                other_train_n = other_train_n(randperm(length(other_train_n)));
                other_train_n = other_train_n(1);        %a randomly selected train number other than the one currently being used for the fit
                curr_trs = find(stim_mat_simple(:, 2) == odor_ni & stim_mat_simple(:, 11) == train_n);
                curr_trs = curr_trs(1:n_trs_to_analyse);
                curr_trs2 = find(stim_mat_simple(:, 2) == odor_ni & stim_mat_simple(:, 11) == other_train_n);
                curr_trs2 = curr_trs2(1:n_trs_to_analyse);                
                ave_dff_resp_mat = mean(dff_data_mat_f(1:(end-5), curr_sig_cells, curr_trs), 3, 'omitnan');
                ave_dff_resp_mat2 = mean(dff_data_mat_f(1:(end-5), curr_sig_cells, curr_trs2), 3, 'omitnan');
                %normalising each dff response trace
                max_vec = max(ave_dff_resp_mat, [], 1);
                max_mat = repmat(max_vec, size(ave_dff_resp_mat, 1), 1);
                ave_dff_resp_mat = ave_dff_resp_mat./max_mat;
                max_vec2 = max(ave_dff_resp_mat2, [], 1);
                max_mat2 = repmat(max_vec2, size(ave_dff_resp_mat2, 1), 1);
                ave_dff_resp_mat2 = ave_dff_resp_mat2./max_mat2;
                
                %Running clustering code
                Z = linkage(ave_dff_resp_mat', 'centroid');
                clust_ids = cluster(Z, 'maxclust', 5);          %grouping cells into a maximum of 5 clusters
                clust_ids_old = clust_ids;
                %giving small clusts the highest clust-id numbers
                n_clusts = max(clust_ids);
                n_cells_list = zeros(n_clusts, 2);
                for clust_n = 1:n_clusts
                    n_cells_list(clust_n, 2) = clust_n;
                    n_cells_list(clust_n, 1) = length(find(clust_ids == clust_n));
                end
                n_cells_list = sortrows(n_cells_list, -1);
                
                for clust_n = 1:n_clusts
                    old_clust_n = n_cells_list(clust_n, 2);
                    curr_cells = find(clust_ids_old == old_clust_n);
                    clust_ids(curr_cells) = clust_n;
                end
                
                
                %re-ordering responses acc to cluster identity
                sorted_resps = [];
                for clust_n = 1:n_clusts
                    curr_cells = find(clust_ids == clust_n);
                    curr_resps = ave_dff_resp_mat(:, curr_cells);       %response vectors for currently clustered cells
                    ave_resp = mean(curr_resps, 2, 'omitnan');
                    
                    corrs = corrcoef([ave_resp, curr_resps]);
                    corrs = corrs(1, 2:end);
                    curr_resps = [corrs; curr_resps];
                    curr_resps = sortrows(curr_resps', -1)';
                    sorted_resps = [sorted_resps, curr_resps(2:end, :)];
                      
                end
                curr_train = stim_mat(curr_trs(1)).rand_trains;
                stim_latency = stim_mat(curr_trs(1)).stim_latency;
                od_pulse_frames = compute_pulse_frames_train(curr_train, frame_time, stim_latency);
                
                figure(1)
                imagesc(sorted_resps', [0, 1])
                colormap(greymap)
                colorbar
                ylabel('cell number')
                set_xlabels_time(1, (frame_time./1000), 3)
                fig_wrapup(1)
                add_stim_bar(1, od_pulse_frames, color_vec(odor_n, :))
                
                %calculating new corrcoeff mat with re-arranged cells
                figure(2)
                corr_mat = corrcoef(sorted_resps);
                imagesc(corr_mat)
                xlabel('cell number')
                ylabel('cell number')
                colorbar
                fig_wrapup(2)
                
                keyboard
                close all
            end
        end
        %keyboard
        
        
        
    end
    
    keyboard
end
       