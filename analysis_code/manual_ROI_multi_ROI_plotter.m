clear all
close all

list_path = 'C:\Data\Code\general_code_old\data_folder_lists\Janelia\13F02Gal4_x_ACh3_choline_sensor.xls';

suppress_plots = 0;
plotting_quant_no_filt = 0;     %1 - only unfiltered traces used for all analysis and plotting - traces included. 0 - filtered traces used for everything.

[del, odor_names1] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
[del, odor_names2] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList_olf2.xls', 1);
odor_names2{3} = 'Butyl acetate';

global color_vec;                
a = colormap('bone');
global greymap
greymap = flipud(a);
fly_n = 0;
script_name = mfilename;

y_ax_traces = 0.8;
y_ax_fit_traces = 0.6;

an_save_path = 'C:\Data\Data\Analysed_data\Analysis_results\PaBaEl_Gamma2\';
force_resave = 1;

n_sec = 2;      %width of moving integration window in s


[del, dir_list] = xlsread(list_path, 1);        %list of Suite2P results directories
n_dirs = size(dir_list, 1);
dataset_list_name = findstr(list_path, 'list_');
dataset_list_name = list_path((dataset_list_name + 5):(end - 4));


%loop to go through all experiment datasets listed in list file
for dir_n = 1:n_dirs
    fly_n = fly_n + 1;

    curr_dir = [dir_list{dir_n, 1}, '\'];
    curr_dir = manage_base_paths(curr_dir, 2)

    tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
    tif_times = tif_times.time_stamps;

    [stim_mat, stim_mat_simple, column_heads, color_vec, good_tr_list, params_orig] = load_params_trains_modular(curr_dir, tif_times);    %reading in trial stimulus parameters after matching time stamps to F traces
    paired_color = color_vec(2, :);
    unpaired_color = color_vec(1, :);
    EL_color = color_vec(3, :);
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

    odor_list_olf2 = unique(stim_mat_simple(:, od_col_ns(2) ) );
    n_odors_olf2 = length(odor_list_olf2);
    odor_dur_list_olf2 = unique(stim_mat_simple(:, dur_col_ns(2) ) );
    n_od_durs_olf2 = length(odor_dur_list_olf2);


    cd(curr_dir);
    tif_name = dir('*.tif');
    if isempty(tif_name) == 1
        continue
    else
    end
    stack_obj = ScanImageTiffReader([curr_dir, tif_name(1).name]);
    [frame_time, zoom, n_chans, PMT_offsets] = SI_tif_info(stack_obj);


    %loading extracted raw fluorescence data matrices written by raw_dff_extractor
    raw_data_mat = load([curr_dir 'extracted_raw_data_mat.mat']);
    raw_data_mat = raw_data_mat.raw_data_mat;           %raw F traces extracted from ROIs
    raw_data_mat_orig = raw_data_mat;
    tif_n_col_n = find_stim_mat_simple_col('matched_tif_n', column_heads);

    %inserting dummy trials (all nans) into raw_data_mat for pairing trials for
    %which no corress .tifs were acquired
    [raw_data_mat, good_tr_list] = match_up_rawmat_matchedtrs(raw_data_mat, stim_mat_simple, tif_n_col_n, good_tr_list);

    n_cells = size(raw_data_mat, 2);

    %calculating dF/F traces from raw data
    filt_time = 0.2;            %in s, the time window for boxcar filter for generating filtered traces
    [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, curr_dir, good_tr_list);
    if plotting_quant_no_filt == 1
        dff_data_mat_f = dff_data_mat;
    else
    end

    del = find(dff_data_mat_f < -1);
    dff_data_mat_f(del) = -1;       %forcing crazy values to sane ones
    
    
    %plotting averaged image with ROIs overlaid
    ave_im_stack = load([curr_dir, 'tr_avg_stack.mat']);
    ave_im_stack = ave_im_stack.ave_stack;
    ROI_mat = load([curr_dir, 'ROI_mat.mat']);
    ROI_mat = ROI_mat.ROI_mat;
    frame = squeeze(ave_im_stack(:, :, 1));
    curr_thresh = median(reshape(frame, 1, []), 'omitnan').*4.5;    %scaling colormap to median pix intensity
    figure(1)
    imagesc(frame, [0, curr_thresh]);
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    colormap('bone')
    figure(2)
    imagesc(frame, [0, curr_thresh]);
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    colormap('bone')
    hold on
    summed_ROI_mat = sum(ROI_mat, 3);
    summed_ROI_mat = imbinarize(summed_ROI_mat, 0.5);
    ROI_mat_sc = summed_ROI_mat.*max(max(frame));
    ROI_obj = imagesc(ROI_mat_sc);
    hold off
    ROI_obj.AlphaData = summed_ROI_mat.*0.5;
    
    
    %plotting traces for each odor stimulus (type)
    
    %replacing nans in stim_mat_simple with 0s.
    stim_mat_simple(isnan(stim_mat_simple)) = 0;
    
    %cascade plot vars
    ypatch_size = 0.5;
    max_y_ax = 0.4;
    x_offset = 0;
    y_offset = 0.7;
    spec_ROI = 1;
    
    %identifying handover trials 
    handover_trs = find(stim_mat_simple(:, dur_col_ns(1)) > 1 & stim_mat_simple(:, dur_col_ns(2)) > 1);
    %looping through olf1 odors, simple trials
    for od_n = 1:length(odor_list_olf1)
        od_ni = odor_list_olf1(od_n);
        curr_trs = find(stim_mat_simple(:, od_col_ns(1)) == od_ni);     %all trials with current odor
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
        stim_frs = stim_frs{1};

        [del1, del2] = intersect(curr_trs, handover_trs);
        if isempty(del2) == 0
            curr_trs(del2) = [];
        else
        end
        
        if isempty(curr_trs) == 1
            continue
        elseif stim_mat_simple(curr_trs(1), dur_col_ns(1)) < 1
            continue
        else 
        end
        keyboard
        fig_h = figure('Name',['olf1 simple od, ', odor_names1{od_ni}]);
        curr_traces = dff_data_mat_f(:, :, curr_trs);
        
        mean_traces = mean(curr_traces, 3, 'omitnan');
        cascade_plot(fig_h, mean_traces, [0.6, 0.6, 0.6], 2, x_offset, y_offset, ypatch_size, max_y_ax);
        set_xlabels_time(fig_h, frame_time, 10)
        fig_wrapup(fig_h, [])
        add_stim_bar(fig_h, stim_frs, color_vec(od_n, :))
        
        %plotting mean +-SE for a single ROI as selected by user
        mean_trace = mean_traces(:, spec_ROI);
        se_trace = std(squeeze(curr_traces(:, spec_ROI, :)), [], 2, 'omitnan')./sqrt(size(curr_traces, 2));
        fig_h = figure('Name',['olf2, simple, ', odor_names2{od_ni}]);
        shadedErrorBar([], mean_trace, se_trace, {'Color', [0.6, 0.6, 0.6]});
        set_xlabels_time(fig_h, frame_time, 10);
        fig_wrapup(fig_h, [])
        add_stim_bar(fig_h, stim_frs, color_vec(od_n, :))
        
    end
    
    %looping through olf2 odors, simple trials
    for od_n = 1:length(odor_list_olf2)
        od_ni = odor_list_olf2(od_n);
        curr_trs = find(stim_mat_simple(:, od_col_ns(2)) == od_ni);     %all trials with current odor
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
        stim_frs = stim_frs{2};

        [del1, del2] = intersect(curr_trs, handover_trs);
        if isempty(del2) == 0
            curr_trs(del2) = [];
        else
        end
        
        if isempty(curr_trs) == 1
            continue
        else
        end
        
        fig_h = figure('Name',['olf2, simple, ', odor_names2{od_ni}]);
        curr_traces = dff_data_mat_f(:, :, curr_trs);
        
        mean_traces = mean(curr_traces, 3, 'omitnan');
        cascade_plot(fig_h, mean_traces, [0.6, 0.6, 0.6], 2, x_offset, y_offset, ypatch_size, max_y_ax);
        set_xlabels_time(fig_h, frame_time, 10);
        fig_wrapup(fig_h, [])
        add_stim_bar(fig_h, stim_frs, color_vec(od_n, :))
        
        %plotting mean +-SE for a single ROI as selected by user
        mean_trace = mean_traces(:, spec_ROI);
        se_trace = std(squeeze(curr_traces(:, spec_ROI, :)), [], 2, 'omitnan')./sqrt(size(curr_traces, 2));
        fig_h = figure('Name',['olf2, simple, ', odor_names2{od_ni}]);
        shadedErrorBar([], mean_trace, se_trace, {'Color', [0.6, 0.6, 0.6]});
        set_xlabels_time(fig_h, frame_time, 10);
        fig_wrapup(fig_h, [])
        add_stim_bar(fig_h, stim_frs, color_vec(od_n, :))
    end
    
    
    %plotting traces for transitioning odor stimuli
    for od_n1 = 1:length(odor_list_olf1)
        od_n1i = odor_list_olf1(od_n1);
        for od_n2 = 1:length(odor_list_olf2)
            od_n2i = odor_list_olf2(od_n2);
            
            curr_trs = find(stim_mat_simple(:, od_col_ns(1)) == od_n1i & stim_mat_simple(:, od_col_ns(2)) == od_n2i );
            curr_trs = intersect(curr_trs, handover_trs);
           
            if isempty(curr_trs) == 1
                continue
            else
            end
            stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
            stim_frs = [stim_frs{1}; stim_frs{2}];
            fig_h = figure('Name',['transitioning, ', odor_names1{od_n1i}, ' - ', odor_names2{od_n2i}]);
            curr_traces = dff_data_mat_f(:, :, curr_trs);
            mean_traces = mean(curr_traces, 3, 'omitnan');
            cascade_plot(fig_h, mean_traces, [0.6, 0.6, 0.6], 2, x_offset, y_offset, ypatch_size, max_y_ax);
            set_xlabels_time(fig_h, frame_time, 10);
            fig_wrapup(fig_h, [])
            add_stim_bar(fig_h, stim_frs, [color_vec(od_n1, :); color_vec(od_n2, :)])
            
            %plotting mean +-SE for a single ROI as selected by user
            mean_trace = mean_traces(:, spec_ROI);
            se_trace = std(squeeze(curr_traces(:, spec_ROI, :)), [], 2, 'omitnan')./sqrt(size(curr_traces, 2));
            fig_h = figure('Name',['olf2, simple, ', odor_names2{od_ni}]);
            shadedErrorBar([], mean_trace, se_trace, {'Color', [0.6, 0.6, 0.6]});
            set_xlabels_time(fig_h, frame_time, 10);
            fig_wrapup(fig_h, [])
            add_stim_bar(fig_h, stim_frs, [color_vec(od_n1, :); color_vec(od_n2, :)])
            
        end
    end
    
    
    keyboard
    close all
end
