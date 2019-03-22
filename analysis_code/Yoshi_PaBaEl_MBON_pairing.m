clear all
close all

dataset_list_paths = [...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_Yoshi_PaBaEl_MBONalpha1_lowUS.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_YoshiPaBaEl_MBON_alpha1_lowUS_backward_ctrl.xls'}...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_Yoshi_PaBaEl_MBON_DAN_gamma1_lowUS_MB085C.xls'}...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_YoshiPaBaEl_MBON_DAN_gamma1_lowUS_MB085C_epoxy.xls'}...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_YoshiPaBaEl_MBON_gamma1_lowUS.xls'}...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_Yoshi_El_THnull_gamma1pedc.xls'}...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_YoshiPaBaEl_MBON_DAN_gamma1pedc_loctite_LED_1ms_0.1Hz.xls'}...     %variable resps, loctite glue, incomplete curing? 4 flies.
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_YoshiPaBaEl_MBON_DAN_gamma1pedc_loctite_LED_1ms.xls'}...           %1 fly, loctite, 0.5 HZ, 1 ms LED pulses.
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\dataset_list_YoshiPaBaEl_MBON_DAN_gamma1pedc_epoxy_LED_1ms_0.1Hz.xls'}...       %6 flies, epoxy, 0.1Hz, 1 ms LED pulses
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\PaBaEl_MBON_DAN_gamma1_lowUS_MB085C_epoxy_short_session.xls'}...                       %4 flies, epoxy, 0.1Hz, 1ms LED pulses, 2 reps pre and post
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\PaBaEl_MBON_DAN_gamma1_lowUS_MB085C_epoxy_short_session_noUS.xls'}... 
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\PaBaEl_MBON_DAN_gamma1_lowUS_MB085C_epoxy_short_session_noUS _shortCS.xls'}... 
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\PaBaEl_MBON_DAN_gamma1_lowUS_MB085C_epoxy_short_session_low_LED.xls'}... 
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\PaBaEl_MBON_DAN_gamma1_lowUS_MB085C_epoxy_short_session_low_LED_1Hz.xls'}... 
                    ];
            
suppress_plots = 0;
[del, odor_names] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
global color_vec;                
a = colormap('bone');
global greymap
greymap = flipud(a);
fly_n = 0;
script_name = mfilename;

for list_n = 1:size(dataset_list_paths, 1)
    curr_dir_list_path = dataset_list_paths{list_n, 1};
    [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
    n_dirs = size(dir_list, 1);
    dataset_list_name = findstr(curr_dir_list_path, 'El_');
    dataset_list_name = curr_dir_list_path((dataset_list_name + 1):(end - 4));
    
    dataset_list_name(1) = [];
    flies_resp_size_mat = [];
    saved_resp_sizes_all = [];
    
    %loop to go through all experiment datasets listed in list file
    for dir_n = 1:n_dirs
       fly_n = fly_n + 1;
              
       saved_an_results.scriptname = mfilename('fullpath');
       curr_dir = [dir_list{dir_n, 1}, '\'];
       curr_dir = manage_base_paths(curr_dir, 2);
       
       tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
       tif_times = tif_times.time_stamps;
       [stim_mat, stim_mat_simple, column_heads, color_vec] = load_params_trains(curr_dir, tif_times);    %reading in trial stimulus parameters after matching time stamps to F traces
       odor_colors = [color_vec(3, :); color_vec(3, :).*0.75; color_vec(2, :)];
       
       %Reading in experimental parameters
        %odor_list = unique(stim_mat_simple(:, 2) );
        odor_list = [3, 10, 11];
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
        [resp_sizes, sig_trace_mat, sig_trace_mat_old, sig_cell_mat, resp_areas] = cal_sig_responses_res(dff_data_mat, stim_mat, stim_mat_simple, curr_dir, frame_time);
        sig_cells = find(sum(sum(sig_cell_mat, 3), 2) > 0);         %list of all cells significant for any odor for any duration
        
        
        %putting in condition that MBON should be significantly responsive
        %to all three odorants in the pre or post trials.
        if sum(sig_cell_mat(:, :, 1), 2) < 3
            disp('warning: not sigresp to all odors')
            %continue
        else
        end
        
        
        figure(1)
        imagesc(squeeze(dff_data_mat_f(:, 1, :))', [0, 1])
        stim_frs = compute_stim_frs(stim_mat, 1, frame_time);
        xlabel('trial n')
        set_xlabels_time(1, frame_time, 10)
        
        fig_wrapup(1, script_name);
        add_stim_bar(1, stim_frs, [0.5, 0.5, 0.5])
        colormap(greymap)
        
        %identifying pre and post pairing trial sets
        pairing_tr = find(stim_mat_simple(:, 13) == 1);
        if isempty(pairing_tr) == 1
            pairing_tr = 15;
        else
        end
        %last_csminus_tr = pairing_tr + length(odor_list) - 1;
        last_csminus_tr = pairing_tr + 1;
        share_path = 'C:\Data\Data\Analysed_data\data_sharing\Glenn_talk_20181121\';
        saved_resp_sizes = [];
        for odor_n = 1:n_odors
            odor_ni = odor_list(odor_n);
            curr_trs = find(stim_mat_simple(:, 2) == odor_ni);
            if odor_n == 1 
                curr_colour = color_vec(3, :);
            elseif odor_n == 2
                curr_colour = color_vec(3, :).*0.75;
            elseif odor_n == 3
                curr_colour = color_vec(2, :);
            end
            
            %plotting pre-trials traces
            figure(2)
            traces_pre = squeeze(dff_data_mat_f(:, 1, curr_trs(curr_trs < pairing_tr)));
            try
                %traces_pre(:, 1:2) = [];                 %not considering first two trials bec of habituation effect
            catch
                keyboard
            end
            
            trace_lengths = size(traces_pre, 1) - sum(isnan(traces_pre(:, 1)));
            trace_lengths = max([trace_lengths, 1]);
            stim_frs = compute_stim_frs(stim_mat, 1, frame_time);
            mean_trace = mean(traces_pre, 2, 'omitnan');
            tr_size = max(mean_trace(stim_frs(1, 1):(stim_frs(1,2))), [], 1);
            fly_resp_size_vec((odor_n - 1).*2 + 1) = tr_size;
            ses = std(traces_pre, [], 2)./sqrt(size(traces_pre, 2));
            %shadedErrorBar([], mean_trace, ses, {'Color', [166, 156, 204]./256}, 1)
            shadedErrorBar([], mean_trace, ses, {':', 'Color', curr_colour, 'lineWidth', 1}, 1);
            axis([0, trace_lengths, -0.25, 1])
            odor_name = odor_names{odor_ni, 1};
            ylabel([odor_name, ' dF/F'])
            hold on           
            
            %plotting post-trials' traces
            traces_post = squeeze(dff_data_mat_f(:, 1, curr_trs(curr_trs > last_csminus_tr)));
            %plot(traces_post, 'Color', [0.6, 0.6, 0.6]);
            stim_frs = compute_stim_frs(stim_mat, 1, frame_time);
            mean_trace = mean(traces_post, 2, 'omitnan');
            tr_size = max(mean_trace(stim_frs(1, 1):(stim_frs(1,2))), [], 1);
            fly_resp_size_vec((odor_n - 1).*2 + 2) = tr_size;
            ses = std(traces_post, [], 2)./sqrt(size(traces_post, 2));
            %shadedErrorBar([], mean_trace, ses, {'Color', [247, 148, 29]./256}, 1)
            shadedErrorBar([], mean_trace, ses, {'Color', curr_colour.*0.75, 'lineWidth', 1}, 1);
            %plot(mean_trace', 'Color', 'k', 'LineWidth', 3);
            axis([0, trace_lengths, -0.25, 2])
            odor_name = odor_names{odor_ni, 1};
            set_xlabels_time(2, frame_time, 10)
            script_name = mfilename;
            fig_wrapup(2, script_name);
            add_stim_bar(2, stim_frs, [0.5, 0.5, 0.5])
            
            %saving data to file for sharing
            metadata.stim_frs = stim_frs;
            metadata.frame_time_s = frame_time;
            
            save([share_path, 'metadata.mat'], 'metadata');
            save([share_path, odor_name, '_pre.mat'], 'traces_pre');
            save([share_path, odor_name, '_post.mat'], 'traces_post');
            
            disp(curr_dir)
            
            %plotting odor response sizes across trials for each fly
            odor_ni = odor_list(odor_n);
            curr_trs = find(stim_mat_simple(:, 2) == odor_ni);
            curr_resps = resp_areas(1, curr_trs);
            
            try
                saved_resp_sizes = pad_n_concatenate(saved_resp_sizes, curr_resps, 1, nan);
            catch 
                keyboard
            end
            
            
            if suppress_plots == 0
                keyboard
            else
            end
            
            
            
            close figure 2
            
        end
        saved_resp_sizes_all = pad_n_concatenate(saved_resp_sizes_all, saved_resp_sizes, 3, nan);
        clear saved_resp_sizes;
        
        flies_resp_size_mat(fly_n, :) = fly_resp_size_vec;
        
        saved_sizes(fly_n).list_name = dataset_list_name;
        saved_sizes(fly_n).resp_sizes = fly_resp_size_vec;
        
        
        
        close figure 1
        
    end
    
    %plotting response sizes across flies for this dataset list
    figure(4)
    marker_colors = [];
    marker_colors(1:2, :) = [repmat(color_vec(3, :), 2, 1)];
    marker_colors(3:4, :) = [repmat(color_vec(3, :).*0.75, 2, 1)];
    marker_colors(5:6, :) = [repmat(color_vec(2, :), 2, 1)];
    
    ylabel('peak dF/F')
    col_pairs = [1, 2; 3, 4; 5, 6];         %this is a list of pairs of columns of points to be connected by lines in the plot
    fig_h4 = scattered_dot_plot(flies_resp_size_mat, 4, 1, 4, 8, marker_colors, 1, col_pairs, [0.75, 0.75, 0.75],...
                            [{'PA pre'}, {'PA post'}, {'BA pre'}, {'BA post'}, {'EL pre'}, {'EL post'}], 1, [0.35, 0.35, 0.35]);
    fig_wrapup(4, script_name)
                        
    %plotting same data as in fig 4 with bars
    line_color = [0.6, 0.6, 0.6];
    pre_color = [78, 54, 255]./256;
    post_color = [255, 130, 0]./256;
    bar_width = 10;
    bar_space = 10;
    
    fig_h5 = bar_line_plot(5, flies_resp_size_mat, line_color, pre_color, post_color, bar_width, bar_space);
    ylabel('peak dF/F')
    xticklabels({'PA', 'BA', 'EL'})
    fig_wrapup(5, script_name)
    
    %statistical testing
    %PA
    [hPA, pPA] = ttest(flies_resp_size_mat(:, 1), flies_resp_size_mat(:, 2))
    %BA
    [hBA, pBA] = ttest(flies_resp_size_mat(:, 3), flies_resp_size_mat(:, 4))
    %EL
    [hEL, pEL] = ttest(flies_resp_size_mat(:, 5), flies_resp_size_mat(:, 6))
    
    %plotting response sizes across repeats for each odor, for all flies.
    figure(7)
    mean_resp_sizes = mean(saved_resp_sizes_all, 3, 'omitnan');
    se_resp_sizes = std(saved_resp_sizes_all, [], 3, 'omitnan')./sqrt(size(saved_resp_sizes_all, 3));
    
    for odor_n = 1:3
        shadedErrorBar([], mean_resp_sizes(odor_n, :), se_resp_sizes(odor_n, :), {'Color', odor_colors(odor_n, :)}, 1);
        hold on
    end
    
    figure(8)
    for odor_n = 1:3
        plot( squeeze(saved_resp_sizes_all(odor_n, :, :)), 'Color', odor_colors(odor_n, :));
        hold on
    end
    
end
%save_path = 'C:\Users\Mehrab\Dropbox (HHMI)\data_sharing\Glenn_talk_2018\slide_30\';
%save([save_path, 'resp_size_mat.mat'], 'flies_resp_size_mat');