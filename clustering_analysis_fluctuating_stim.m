clear all
close all

direc_lists_mat =  [{'C:\Data\Data\Analysed_data\dataset_list_fluc_stim_somas_20171226.xls'}...
                    
                   ]; 

save_path = 'C:\Data\Analysis_results\train_clustering\';
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

[del, odor_names] = xlsread('C:\Data\Data\Analysed_data\odor_names_20161108.xls', 1);




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
                ave_dff_resp_mat = mean(dff_data_mat_f(:, curr_sig_cells, curr_trs), 3, 'omitnan');
                ave_dff_resp_mat2 = mean(dff_data_mat_f(:, curr_sig_cells, curr_trs2), 3, 'omitnan');
                %normalising each dff response trace
                max_vec = max(ave_dff_resp_mat, [], 1);
                max_mat = repmat(max_vec, size(ave_dff_resp_mat, 1), 1);
                ave_dff_resp_mat = ave_dff_resp_mat./max_mat;
                max_vec2 = max(ave_dff_resp_mat2, [], 1);
                max_mat2 = repmat(max_vec2, size(ave_dff_resp_mat2, 1), 1);
                ave_dff_resp_mat2 = ave_dff_resp_mat2./max_mat2;
                
                
                %Running clustering algorithm
                X = ave_dff_resp_mat';       %activity data matrix - cells x frames, with trials concatenated

                thr = 0.8;              %threshold of kmeans++ runs a pair of cells must co-occur in to be considered for clustering
                no_iters = 500;         %no of times meta_k_means_tank calls k_means_plusplus; set to 100 for preliminary corrcut scan
                clustat = [];
                un_clust = [];

                %loop to try out various values of corrcut
                ct_range = 0.65:0.05:0.9;
                for corrcut = 0.65:0.05:0.9
                    %corrcut                             %corrcut is the threshold inter-cluster correlation coefficient above which they are merged
                    [a b c] = meta_k_means_tank1(X, 'correlation', corrcut, thr, no_iters);


                    clustcorr = [];
                    clustnum = [];
                    for i = 1 : size(a, 1)
                        temp = a{i,1};
                        clustnum = [clustnum length(temp)];
                        tempcorr = 0;
                        count = 1;
                        for x = 1 : length(temp)-1
                            for y = x + 1 : length(temp)
                                tempcorr = tempcorr + ccorr(temp(x), temp(y));
                                count = count + 1;
                            end
                        end
                        clustcorr = [clustcorr, tempcorr/count];
                    end
                    clustat = [clustat; [mean(clustnum) (sum(clustcorr.*clustnum)/sum(clustnum))]];
                end

                %identifying best value of corrcut
                x = clustat(:, 1);
                x = (x - min(x))/(max(x)-min(x));
                y = clustat(:, 2);
                y = (y - min(y))/(max(y)-min(y));
                [null, corrcuti] = max(x.*y);

                corrcut = ct_range(corrcuti);

                clear null

                %running meta-k-means to identify final clusters with
                %well-chosen value of corrcut (threshold corrcoef to fuse meta-clusters)
                no_iters = 1000;            
                [a, b, c] = meta_k_means_tank1(X, 'correlation', corrcut, thr, no_iters);

                %saving results of clustering to file 
                clust_results = {a, b, c};
                mkdir([save_path, dir_list_name])
                save([save_path, dir_list_name, '.mat'], 'clust_results');
                
                beep

                %building list of cells and their group numbers
                no_clusts = size(a, 1);
                cell_gp_vec = zeros(n_cells, 1);
                for c_num = 1:no_clusts
                    c_list = a{c_num, 1};
                    cell_gp_vec(c_list, 1) = c_num;
                end

                clear no_clusts
                clear c_num
                clear c_list

                cell_gp_vec = [cell_gp_vec, (1:1:n_cells)'];
                cell_gp_vec_orig = cell_gp_vec;

                %sorting cells within groups as per their corrcoeffs with
                %group-averaged trace
                c_lengths = [];
                c_list_f = cell_gp_vec_orig(:, 1);                  %full list of all cells' cluster numbers
                
                response_matrix_o = [];
                response_matrix_o_20 = [];
                %accounting for un-clustered cells
                c_list = find(c_list_f == 0);
                un_clust = [un_clust; length(c_list), n_cells];
                for c_no = 1:length(c_list)
                    c_noi = c_list(c_no);

                end


                %sorting response matrix by corrcoef with cluster mean
                for clust_no = 1:max(cell_gp_vec_orig(:, 1))
                    c_list = find(c_list_f == clust_no);            %cell numbers in old list that belong to current cluster

                    %condition to skip clusters with < 5 cells in them
%                     if length(c_list) < 5
%                         for c_no = 1:length(c_list)
%                             c_noi = c_list(c_no);
%                         end
% 
% 
%                         continue
%                     else
%                     end                
                    c_lengths = [c_lengths; length(c_list)]; 
                    cr_list = c(c_list, clust_no);                  %list of corrcoeffs of cells in c_list, with their own cluster's averaged trace 
                    curr_trace_mat = saved_traces_flies(:, c_list, long_dur_n)';       %matrix of activity data for each cell that belongs to this cluster
                    curr_trace_mat_20 = saved_traces_flies(:, c_list, (long_dur_n-1) )';             %matrix of activity data for 20s odor trials 

                    curr_trace_mat = [cr_list, curr_trace_mat];     %concatenating corrcoefs to activity traces to sort both together
                    curr_trace_mat_20 = [cr_list, curr_trace_mat_20];

                    curr_trace_mat = sortrows(curr_trace_mat, -1);      %sorting rows by first column ie corrcoefs
                    curr_trace_mat_20 = sortrows(curr_trace_mat_20, -1);

                    response_matrix_o = [response_matrix_o; curr_trace_mat(:, 2:end)];
                    response_matrix_o_20 = [response_matrix_o_20; curr_trace_mat_20(:, 2:end)];
                end

                max_mat = repmat(nanmax(response_matrix_o, [], 2), 1, size(response_matrix_o, 2));
                response_matrix_o_norm = response_matrix_o./max_mat;

                max_mat = repmat(nanmax(response_matrix_o_20, [], 2), 1, size(response_matrix_o_20, 2));
                response_matrix_o_norm_20 = response_matrix_o_20./max_mat;

                if isempty(response_matrix_o_norm) == 1
                    continue
                else
                end
                response_matrix_o_norm_s = tsmovavg_m(response_matrix_o_norm, 's', 10, 2);
                response_matrix_o_norm_s(:, 1:9) = response_matrix_o_norm(:, 1:9);

                response_matrix_o_norm_20_s = tsmovavg_m(response_matrix_o_norm_20, 's', 10, 2);
                response_matrix_o_norm_20_s(:, 1:9) = response_matrix_o_norm_20(:, 1:9);

                %calculating new corrcoeff mat with re-arranged cells
                curr_color = color_vec(2, :);
                fig_h = figure(1);
                imagesc(response_matrix_o_norm_s, [0, 1]);
                stim_frs_saved = [stim_frame, stim_end_fr];
                colormap(greymap)
                set_xlabels_time(1, frame_time, .5)
                ylabel('sig. cell-odor pairs')
                xlabel('time (s)')
                set(fig_h, 'Position', [100, 100, 100 + plot_width, 100 + plot_height]);
                fig_wrapup(1);
                add_stim_bar(1, stim_frs_saved, curr_color)





                
            end
        end
        %keyboard
        
        
        
    end
    
    keyboard
end
       