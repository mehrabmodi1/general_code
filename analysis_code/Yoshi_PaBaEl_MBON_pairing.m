clear all
close all

dataset_list_paths = [...
                    {'C:\Data\Code\general_code\data_folder_lists\Janelia\dataset_list_Yoshi_PaBaEl_MBONalpha1_lowUS.xls'};...


                ];
            
suppress_plots = 1;
[del, odor_names] = xlsread('C:\Data\Code\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
global color_vec;                
color_vec = load('C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Code\std_color_vec.txt');
a = colormap('bone');
global greymap
greymap = flipud(a);
fly_n = 0;

for list_n = 1:size(dataset_list_paths, 1)
    curr_dir_list_path = dataset_list_paths{list_n, 1};
    [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
    n_dirs = size(dir_list, 1);
    dataset_list_name = findstr(curr_dir_list_path, 'El_');
    dataset_list_name = curr_dir_list_path((dataset_list_name + 1):(end - 4));
    dataset_list_name(1) = [];
    flies_resp_size_mat = [];
    
    %loop to go through all experiment datasets listed in list file
    for dir_n = 1:n_dirs
       fly_n = fly_n + 1;
              
       saved_an_results.scriptname = mfilename('fullpath');
       curr_dir = [dir_list{dir_n, 1}, '\'];
       tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
       tif_times = tif_times.time_stamps;
       [stim_mat, stim_mat_simple, column_heads, color_vec] = load_params_trains(curr_dir, tif_times);    %reading in trial stimulus parameters after matching time stamps to F traces
       
       %Reading in experimental parameters
        odor_list = unique(stim_mat_simple(:, 2) );
        n_odors = length(odor_list);
        odor_dur_list = unique(stim_mat_simple(:, 3) );
        n_od_durs = length(odor_dur_list);
        n_trains = max(stim_mat_simple(:, 11));
        saved_an_results.odor_list = odor_list;
        saved_an_results.odor_dur_list = odor_dur_list;
        fly_resp_size_vec = zeros(1, length(odor_list).*2) + nan;
        
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
        del = find(dff_data_mat_f < -1);
        dff_data_mat_f(del) = -1;       %forcing crazy values to sane ones
        
        %identifying significantly responsive cells
        [resp_sizes, sig_trace_mat, sig_trace_mat_old, sig_cell_mat] = cal_sig_responses_res(dff_data_mat, stim_mat, stim_mat_simple, curr_dir, frame_time);
        sig_cells = find(sum(sum(sig_cell_mat, 3), 2) > 0);         %list of all cells significant for any odor for any duration
        
        figure(1)
        imagesc(squeeze(dff_data_mat_f(:, 1, :))', [0, 1])
        stim_frs = compute_stim_frs(stim_mat, 1, frame_time);
        xlabel('trial n')
        set_xlabels_time(1, frame_time, 10)
        fig_wrapup(1);
        add_stim_bar(1, stim_frs, [0.5, 0.5, 0.5])
        colormap(greymap)
        
        for odor_n = 1:n_odors
            odor_ni = odor_list(odor_n);
            curr_trs = find(stim_mat_simple(:, 2) == odor_ni);
                        
            %plotting pre-trials traces
            figure(2)
            plot(squeeze(dff_data_mat_f(:, 1, curr_trs(curr_trs<16))), 'Color', [0.6, 0.6, 0.6]);
            hold on
            stim_frs = compute_stim_frs(stim_mat, 1, frame_time);
            mean_trace = mean(squeeze(dff_data_mat_f(:, 1, curr_trs(curr_trs<15))), 2);
            tr_size = mean(mean_trace(stim_frs(1, 1):stim_frs(1,2)), 1);
            fly_resp_size_vec((odor_n - 1).*2 + 1) = tr_size;
            plot(mean_trace', 'Color', 'k', 'LineWidth', 3);
            axis([0, 800, -0.25, 1])
            odor_name = odor_names{odor_ni, 1};
            ylabel([odor_name, ' pre-trs (dF/F)'])
            set_xlabels_time(2, frame_time, 10)
            fig_wrapup(2);
            add_stim_bar(2, stim_frs, [0.5, 0.5, 0.5])
            
            
            %plotting post-trials traces
            figure(3)
            plot(squeeze(dff_data_mat_f(:, 1, curr_trs(curr_trs > 18))), 'Color', [0.6, 0.6, 0.6]);
            hold on
            stim_frs = compute_stim_frs(stim_mat, 1, frame_time);
            mean_trace = mean(squeeze(dff_data_mat_f(:, 1, curr_trs(curr_trs > 18))), 2);
            tr_size = mean(mean_trace(stim_frs(1, 1):stim_frs(1,2)), 1);
            fly_resp_size_vec((odor_n - 1).*2 + 2) = tr_size;
            plot(mean_trace', 'Color', 'k', 'LineWidth', 3);
            axis([0, 800, -0.25, 1])
            odor_name = odor_names{odor_ni, 1};
            ylabel([odor_name, ' post-trs (dF/F)'])
            set_xlabels_time(3, frame_time, 10)
            fig_wrapup(3);
            add_stim_bar(3, stim_frs, [0.5, 0.5, 0.5])
            
            %del = input('press enter for next odor.');
            
            
            close figure 2
            close figure 3
            
            
            
        end
        
        flies_resp_size_mat(fly_n, :) = fly_resp_size_vec;
        
        saved_sizes(fly_n).list_name = dataset_list_name;
        saved_sizes(fly_n).resp_sizes = fly_resp_size_vec;
        %keyboard
        
        
        close figure 1
        
    end
    
    %plotting response sizes across flies for this dataset list
    figure(4)
    marker_colors = [];
    for odor_n = 1:n_odors
        if odor_n == 1
            marker_colors = [marker_colors; repmat(color_vec(2, :).*0.75, 2, 1)];
        else 
            marker_colors = [marker_colors; repmat(color_vec(odor_n, :), 2, 1)];
       
        end
    
    end
    
    fig_h4 = scattered_dot_plot(flies_resp_size_mat, 4, 2, 4, 8, marker_colors, 1, 0, [0.75, 0.75, 0.75], [{'PA pre'}, {'PA post'}, {'BA pre'}, {'BA post'}, {'EL pre'}, {'EL post'}], 1, [0.65, 0.65, 0.65]);
    
    %statistical testing
    %PA
    [hPA, pPA] = ttest(flies_resp_size_mat(:, 1), flies_resp_size_mat(:, 2))
    %BA
    [hBA, pBA] = ttest(flies_resp_size_mat(:, 3), flies_resp_size_mat(:, 4))
    %EL
    [hEL, pEL] = ttest(flies_resp_size_mat(:, 5), flies_resp_size_mat(:, 6))
    keyboard 
    
end
        