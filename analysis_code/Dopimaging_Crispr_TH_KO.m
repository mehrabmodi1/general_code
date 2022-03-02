clear all
close all

dataset_list_paths = [...
                      
                      {'C:\Data\Code\general_code\data_folder_lists\Janelia\13F02LexADANsensor_Cas9_THgRNAinDANG2Ap1.xls'};...
                      
                      ];
            
suppress_plots = 1;
cell_n = 1;

global color_vec;                
a = colormap('bone');
global greymap
greymap = flipud(a);
fly_n = 0;
script_name = mfilename;

y_ax_traces = 0.8;
y_ax_fit_traces = 0.6;

an_save_path = 'C:\Data\Data\Analysed_data\Analysis_results\Crispr_KOs\DA_img_20211230\';
force_resave = 1;

integ_win_s = 8;        %width of pulse integration window for response quantification
n_trs = 2;

plotting_quant_no_filt = 0;    %manual toggle to decide whether to use filtered traces for plots, quantification

for list_n = 1:size(dataset_list_paths, 1)
    
    curr_dir_list_path = dataset_list_paths{list_n, 1};
    [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
    n_dirs = size(dir_list, 1);
    dataset_list_name = findstr(curr_dir_list_path, 'list_');
    dataset_list_name = curr_dir_list_path((dataset_list_name + 5):(end - 4));
    
    
    set_list_type = 0;  %simple stimulus case
    
    
    saved_traces_all = [];
    saved_odor_resps_all = [];
    %loop to go through all experiment datasets listed in list file
    for dir_n = 1:n_dirs
        fly_n = fly_n + 1;
             
        curr_dir = [dir_list{dir_n, 1}, '\'];
        curr_dir = manage_base_paths(curr_dir, 2)
       
        tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
        tif_times = tif_times.time_stamps;
        
        cd(curr_dir);
        tif_name = dir('*.tif');
        if isempty(tif_name) == 1
            continue
        else
        end
        stack_obj = ScanImageTiffReader([curr_dir, tif_name(1).name]);
        [frame_time, zoom, n_chans, PMT_offsets] = SI_tif_info(stack_obj);
        
        %Plotting ave image with ROIs overlaid
        stack = stack_obj.data();
        stack = permute(stack,[2 1 3]);
        ave_fr = mean(stack, 3).^2;
        
        ROI_mat = load([curr_dir, 'ROI_mat.mat']);
        %computing ROI_mat centroids according to ROI order used for data extraction
        ROI_mat_labels = ROI_mat.ROI_mat;
        for ROI_n = 1:size(ROI_mat_labels, 3)
            curr_ROI = squeeze(ROI_mat_labels(:, :, ROI_n));
            curr_centroid = regionprops(curr_ROI, 'centroid');
            ROI_centroids(ROI_n, :) = curr_centroid.Centroid;
        end
        
        %plotting labelled ROIs       
        ROI_mat_fr = sum(ROI_mat.ROI_mat, 3);
        se = strel('disk', 5);
        ROI_mat_fr_sub = im2bw(imerode(ROI_mat_fr, se), 0.5);
        ROI_mat_fr = ROI_mat_fr - ROI_mat_fr_sub;       %computing ROI outlines
        figure(1)
        plot_frame(ave_fr, 3, ROI_mat_fr, ROI_centroids);
        
        if exist([an_save_path, 'correct_ROI_ordering_fly_', num2str(dir_n), '.mat']) ~= 2
            correct_ROI_ordering = input('manually specify correct ordering of ROIs');
            save([an_save_path, 'correct_ROI_ordering_fly_', num2str(dir_n), '.mat'], 'correct_ROI_ordering');
        else
            correct_ROI_ordering = load([an_save_path, 'correct_ROI_ordering_fly_', num2str(dir_n), '.mat']);
            correct_ROI_ordering = correct_ROI_ordering.correct_ROI_ordering;
        end
        
        %reordering and re-plotting ROI image
        ROI_centroids_unord = ROI_centroids;
        ROI_centroids = zeros(size(ROI_centroids_unord, 1), 2);
        for ROI_n = 1:size(ROI_centroids, 1)
            ROI_n_old = correct_ROI_ordering(ROI_n);
            ROI_centroids(ROI_n, :) = ROI_centroids_unord(ROI_n_old, :); 
        end
        figure(1)
        plot_frame(ave_fr, 3, ROI_mat_fr, ROI_centroids);
        %keyboard
        
        
        [stim_mat, stim_mat_simple, column_heads, color_vec, good_tr_list, params_orig, PID_traces] = load_params_trains_modular(curr_dir, tif_times, frame_time);    %reading in trial stimulus parameters after matching time stamps to F traces
        
        odor_names1 = stim_mat.odourNames;
        odor_names2 = stim_mat.odourNames_olf2;
        PA_odn_olf1 = od_name_lookup(odor_names1, 'Pentyl acetate');
        PA_odn_olf2 = od_name_lookup(odor_names2, 'Pentyl acetate');
        
        PA_color = color_vec(2, :);
        EL_color = color_vec(3, :);
        
        PA_color_KO = PA_color.*0.6;
        EL_color_KO = EL_color.*0.6;
        %mean_color = [0.5, 0.83, 0.98];
        %mean_color = [149, 200, 216]./256;
        mean_color = ([0, 49, 152]./256).*1.5;
        stim_mat_simple_nonans = stim_mat_simple;
        stim_mat_simple_nonans(isnan(stim_mat_simple)) = 0;
                
        %identifying stim_mat_simple col numbers
        led_on_col_n = find_stim_mat_simple_col('led_on', column_heads);            %identifying relevant column number in stim_mat_simple
        od_olf1_col_n = find_stim_mat_simple_col('odor_n', column_heads);           %identifying relevant column number in stim_mat_simple
        od_olf2_col_n = find_stim_mat_simple_col('odour_olf2', column_heads);       %identifying relevant column number in stim_mat_simple
        dur_olf1_col_n = find_stim_mat_simple_col('duration', column_heads);        %identifying relevant column number in stim_mat_simple
        dur_olf2_col_n = find_stim_mat_simple_col('duration_olf2', column_heads);   %identifying relevant column number in stim_mat_simple
        od_col_ns = [od_olf1_col_n, od_olf2_col_n];
        dur_col_ns = [dur_olf1_col_n, dur_olf2_col_n];
        od_durs = unique(stim_mat_simple(:, dur_col_ns(2)));
        od_durs(isnan(od_durs)) = [];
        odn_list_olf2 = unique(stim_mat_simple(:, od_col_ns(2)));
        odn_list_olf2(isnan(odn_list_olf2)) = [];
                
        %Reading in experimental parameters
        odor_list_olf1 = unique(stim_mat_simple(:, od_col_ns(1) ) );
        n_odors_olf1 = length(odor_list_olf1);
        odor_dur_list_olf1 = unique(stim_mat_simple(:, dur_col_ns(1) ) );
        n_od_durs_olf1 = length(odor_dur_list_olf1);
         
        integ_win = round(integ_win_s./frame_time); %integration window for response quantification in frames
       
        %loading extracted raw fluorescence data matrices written by raw_dff_extractor
        raw_data_mat = load([curr_dir 'extracted_raw_data_mat.mat']);
        raw_data_mat = raw_data_mat.raw_data_mat;           %raw F traces extracted from ROIs
        raw_data_mat_orig = raw_data_mat;
        tif_n_col_n = find_stim_mat_simple_col('matched_tif_n', column_heads);
        n_cells = size(raw_data_mat, 2);
    
        n_cells = size(raw_data_mat, 2);
        
        %re-organizing ROI order according to manually specified correct ordering
        raw_data_mat_unord = raw_data_mat;
        raw_data_mat = zeros(size(raw_data_mat, 1), size(raw_data_mat, 2), size(raw_data_mat, 3));
        for ROI_n = 1:size(ROI_centroids, 1)
            ROI_n_old = correct_ROI_ordering(ROI_n);
            raw_data_mat(:, ROI_n, :) = raw_data_mat_unord(:, ROI_n_old, :); 
        end
        
        %calculating dF/F traces from raw data
        filt_time = 0.2;            %in s, the time window for boxcar filter for generating filtered traces
        
        [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, curr_dir, good_tr_list);
        
        if plotting_quant_no_filt == 1
            dff_data_mat_f = dff_data_mat;
        else
        end
          
        if dir_n == 1
            dff_data_mat_f(:, 6, :) = [];
        else
        end
         
        del = find(dff_data_mat_f < -1);
        dff_data_mat_f(del) = -1;       %forcing crazy values to sane ones
       
        %replacing nans in stim_mat_simple with 0s to make lookup easier
        stim_mat_simple_orig = stim_mat_simple;
        stim_mat_simple(isnan(stim_mat_simple)) = 0;
        
        %computing and logging mean odor resp trace and resp size for each odor, for each ROI 
        saved_od_resps = [];
        saved_traces = [];
        for odor_n = 1:length(odor_list_olf1)
            odor_ni = odor_list_olf1(odor_n);
            curr_trs = find(stim_mat_simple(:, od_col_ns(1)) == odor_ni);
            
            %computing and logging response sizes from unfiltered data
            stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
            integ_frs = [stim_frs{1}(1), (stim_frs{1}(1) + round(integ_win_s./frame_time))];
            curr_resps = mean(dff_data_mat(integ_frs(1):integ_frs(2), :, curr_trs(1:n_trs)), 1, 'omitnan');      %response sizes over odor delivery window for each trial
            mean_resps = squeeze(mean(curr_resps, 3, 'omitnan'));    %resps for each ROI, averaged across trials
            saved_od_resps = pad_n_concatenate(saved_od_resps, mean_resps', 2, nan);  %collating responses for all ROIs across odors
            
           
            
            %plotting dF/F response image
            highlight_resp_pix(5, stack, integ_frs, 1, frame_time, 0.18, 0);
                        
            %computing and logging mean response trace for all ROIs
            mean_traces = squeeze(mean(dff_data_mat_f(:, :, curr_trs(1:n_trs)), 3, 'omitnan'));
            saved_traces = pad_n_concatenate(saved_traces, mean_traces, 3, nan);
            
        end
        saved_odor_resps_all = pad_n_concatenate(saved_odor_resps_all, saved_od_resps, 3, nan);
        saved_traces_all = pad_n_concatenate(saved_traces_all, saved_traces, 4, nan);
        
        stim_mat(1).odours
        keyboard
    end
    
    
    %plotting mean responses for KO and ctrl datasets
    figure(2)
    set(gcf, 'Name', 'reponses to PA');
    resp_mat_KO = [squeeze(saved_odor_resps_all(:, 1, 1:4))];
    resp_mat_ctrl = [squeeze(saved_odor_resps_all(:, 1, 5:8))];
    resp_mat = [resp_mat_ctrl(1, :)', resp_mat_KO(1, :)', resp_mat_ctrl(2, :)', resp_mat_KO(2, :)', resp_mat_ctrl(3, :)', resp_mat_KO(3, :)', resp_mat_ctrl(4, :)', resp_mat_KO(4, :)'...
                    resp_mat_ctrl(5, :)', resp_mat_KO(5, :)', resp_mat_ctrl(6, :)', resp_mat_KO(6, :)', resp_mat_ctrl(7, :)', resp_mat_KO(7, :)'];
    markercolor = [0.75, 0.75, 0.75];
    %xlabels = [{'nogRNA'}, {'THgRNA'}, {'nogRNA'}, {'THgRNA'}, {'nogRNA'}, {'THgRNA'}, {'nogRNA'}, {'THgRNA'}, {'nogRNA'}, {'THgRNA'}, {'nogRNA'}, {'THgRNA'}, {'nogRNA'}, {'THgRNA'}];
    xlabels = [{'1ct'}, {'1KO'}, {'2ct'}, {'2KO'}, {'3ct'}, {'3KO'}, {'4ct'}, {'4KO'}, {'5ct'}, {'5KO'}, {'6ct'}, {'6KO'}, {'7ct'}, {'7KO'},];
    fig_h = scattered_dot_plot_ttest(resp_mat, 2, .6, 2, 4, markercolor, 1, [], [], xlabels, 2, [0, 0, 0], 2, 0.05, 0, 1, 'force_mean');
    
    %plotting mean responses for KO and ctrl datasets
    figure(3)
    set(gcf, 'Name', 'reponses to EL');
    resp_mat_KO = [squeeze(saved_odor_resps_all(:, 2, 1:4))];
    resp_mat_ctrl = [squeeze(saved_odor_resps_all(:, 2, 5:8))];
    resp_mat = [resp_mat_ctrl(1, :)', resp_mat_KO(1, :)', resp_mat_ctrl(2, :)', resp_mat_KO(2, :)', resp_mat_ctrl(3, :)', resp_mat_KO(3, :)', resp_mat_ctrl(4, :)', resp_mat_KO(4, :)'...
                    resp_mat_ctrl(5, :)', resp_mat_KO(5, :)', resp_mat_ctrl(6, :)', resp_mat_KO(6, :)', resp_mat_ctrl(7, :)', resp_mat_KO(7, :)'];
    markercolor = [0.75, 0.75, 0.75];
    %xlabels = [{'nogRNA'}, {'THgRNA'}, {'nogRNA'}, {'THgRNA'}, {'nogRNA'}, {'THgRNA'}, {'nogRNA'}, {'THgRNA'}, {'nogRNA'}, {'THgRNA'}, {'nogRNA'}, {'THgRNA'}, {'nogRNA'}, {'THgRNA'}];
    xlabels = [{'1ct'}, {'1KO'}, {'2ct'}, {'2KO'}, {'3ct'}, {'3KO'}, {'4ct'}, {'4KO'}, {'5ct'}, {'5KO'}, {'6ct'}, {'6KO'}, {'7ct'}, {'7KO'},];
    fig_h = scattered_dot_plot_ttest(resp_mat, 3, .6, 2, 4, markercolor, 1, [], [], xlabels, 2, [0, 0, 0], 2, 0.05, 0, 1, 'force_mean');
    
    
    %plotting mean response traces for all compartments
    stim_frs = stim_frs{1};
    pooled_traces = [];
    for cpt_n = 1:size(ROI_mat.ROI_mat, 3)
        h = figure();
        set(gcf, 'Name', ['cpt ', num2str(cpt_n), ' responses.']);            
        for od_n = 1:size(odor_list_olf1, 1)
            curr_od_name = odor_names1{odor_list_olf1(od_n)};
            subplot(1, size(odor_list_olf1, 1), od_n);
            
            %plotting ctrl trace
            curr_traces = squeeze(saved_traces_all(:, cpt_n, od_n, 5:8));
            mean_trace = mean(curr_traces, 2, 'omitnan');
            se_trace = std(curr_traces, [], 2, 'omitnan')./sqrt(size(curr_traces, 2));
            shadedErrorBar([], mean_trace, se_trace, {'Color', [0.7, 0.7, 0.7]});
            hold on
            
            %plotting KO trace
            curr_traces = squeeze(saved_traces_all(:, cpt_n, od_n, 1:4));
            mean_trace = mean(curr_traces, 2, 'omitnan');
            se_trace = std(curr_traces, [], 2, 'omitnan')./sqrt(size(curr_traces, 2));
            shadedErrorBar([], mean_trace, se_trace, {'Color', 'k'});
            hold off
            
            set_xlabels_time(h, frame_time, 10)
            fig_wrapup(h, [], [75, 180], 0.6)
            add_stim_bar(h, stim_frs, color_vec(od_n, :));
            
        end
        
        pooled_traces = pad_n_concatenate(pooled_traces, mean_trace, 2, nan);
        
    end
    
    h = figure();
    set(gcf, 'Name', ['all cpt mean responses.']);     
    cascade_plot(h, pooled_traces, [0.7, 0.7, 0.7], 2, 0, 2, 1, .4, 0, 0)
    set_xlabels_time(h, frame_time, 10)
    fig_wrapup(h, [], [250, 50], 0.6)
    add_stim_bar(h, stim_frs, [0, 0, 0]);
    keyboard
   
    
end


function [frame_obj, ROI_obj] = plot_frame(frame, curr_thresh, ROI_mat, ROI_centroids)
    
    if median(reshape(frame, 1, []), 'omitnan') ~= 0
        frame_obj = imagesc(frame, [0, curr_thresh.*median(reshape(frame, 1, []), 'omitnan')]);
    elseif median(reshape(frame, 1, []), 'omitnan') == 0
        frame_obj = imagesc(frame, [0, 1]);
    end
    
    
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    colormap('gray')
    hold on
    ROI_mat_sc = ROI_mat.*max(max(frame));
    ROI_obj = imagesc(ROI_mat_sc);
    ROI_obj.AlphaData = ROI_mat.*0.5;
    
    %adding centroid labels
    for ROI_n = 1:size(ROI_centroids, 1)
        curr_cent = ROI_centroids(ROI_n, :);
        text(curr_cent(1), curr_cent(2), num2str(ROI_n));
    end
    hold off
    
   
end