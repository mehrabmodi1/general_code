clear all
close all

dataset_list_paths = [ ...
                         {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\d5HT1b_similar_od_handovers.xls'};...
                      
                      ];
            
suppress_plots = [];
[del, odor_names1] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
[del, odor_names2] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList_olf2.xls', 1);
odor_names2{3} = 'Butyl acetate';
od_names = [{'PA'}, {'BA'}, {'EL'}];
            

global color_vec;                
a = colormap('bone');
global greymap
greymap = flipud(a);
fly_n = 0;
script_name = mfilename;

y_ax_traces = 0.8;
y_ax_fit_traces = 0.6;
saved_long_traces = 0;
all_sig_frs = [];
for list_n = 1:size(dataset_list_paths, 1)
    curr_dir_list_path = dataset_list_paths{list_n, 1};
    [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
    n_dirs = size(dir_list, 1);
    dataset_list_name = findstr(curr_dir_list_path, 'list_');
    dataset_list_name = curr_dir_list_path((dataset_list_name + 5):(end - 4));
    
    
    %loop to go through all experiment datasets listed in list file
    all_resp_traces = [];
    for dir_n = 1:n_dirs
        fly_n = fly_n + 1;
              
        saved_an_results.scriptname = mfilename('fullpath');
        curr_dir = [dir_list{dir_n, 1}, '\'];
        curr_dir = manage_base_paths(curr_dir, 3);
        curr_dir = [curr_dir, '\1\'];
       
        tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
        tif_times = tif_times.time_stamps;
        cd(curr_dir);
        tif_name = dir('*.tif');
        stack_obj = ScanImageTiffReader([curr_dir, tif_name(1).name]);
        [frame_time, zoom, n_chans, PMT_offsets] = SI_tif_info(stack_obj);
        
        [stim_mat, stim_mat_simple, column_heads, color_vec, good_tr_list, params_orig] = load_params_trains_modular(curr_dir, tif_times, frame_time);    %reading in trial stimulus parameters after matching time stamps to F traces
        
       
        %Reading in experimental parameters
        odor_list_olf1 = unique(stim_mat_simple(:, 1) );
        n_odors_olf1 = length(odor_list_olf1);
        odor_dur_list_olf1 = unique(stim_mat_simple(:, 2) );
        n_od_durs_olf1 = length(odor_dur_list_olf1);
        n_trains_olf1 = max(stim_mat_simple(:, 8));
        
        odor_list_olf2 = unique(stim_mat_simple(:, 9) );
        n_odors_olf2 = length(odor_list_olf2);
        odor_dur_list_olf2 = unique(stim_mat_simple(:, 10) );
        n_od_durs_olf2 = length(odor_dur_list_olf2);
        n_trains_olf2 = max(stim_mat_simple(:, 17));
        
        
        
        
        
        
        %loading extracted raw fluorescence data matrices written by raw_dff_extractor
        raw_data_mat = load([curr_dir 'extracted_raw_data_mat.mat']);
        raw_data_mat = raw_data_mat.raw_data_mat;           %raw F traces extracted from ROIs
        raw_data_mat_orig = raw_data_mat;
        tif_n_col_n = find_stim_mat_simple_col('matched_tif_n', column_heads);
        
        %inserting dummy trials (all nans) into raw_data_mat for pairing trials for
        %which no corress .tifs were acquired
        raw_data_mat = match_up_rawmat_matchedtrs(raw_data_mat, stim_mat_simple, tif_n_col_n);
        
        %dumping data from manually identified, z-drifted trials
        bad_tr_list = 1:size(raw_data_mat, 3);
        bad_tr_list(good_tr_list) = [];
        raw_data_mat(:, :, bad_tr_list) = nan;
        
        n_cells = size(raw_data_mat, 2);
        
        %calculating dF/F traces from raw data
        filt_time = .5;            %in s, the time window for boxcar filter for generating filtered traces
        [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, curr_dir, []);
        
        
        del = find(dff_data_mat_f < -1);
        dff_data_mat_f(del) = -1;       %forcing crazy values to sane ones
        
       
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
        odn_list_olf1 = unique(stim_mat_simple(:, od_col_ns(1)));
        
        %replacing Nans in stim_mat_simple with zeroes
        stim_mat_simple_orig = stim_mat_simple;
        del = isnan(stim_mat_simple);
        stim_mat_simple(del) = 0;
        
        %replacing small vals in olf1_dur column with zeroes
        del = find(stim_mat_simple(:, dur_col_ns(1)) < 1);
        stim_mat_simple(del, dur_col_ns(1)) = 0;
       
        [resp_sizes, sig_trace_mat, sig_cell_mat, sig_cell_mat_key, resp_areaundercurves] = cal_sig_responses_res_modular(dff_data_mat, stim_mat, stim_mat_simple, frame_time, od_col_ns, dur_col_ns);
        
        union_sig_cell_vec = sum(sig_cell_mat, 2);
        union_sig_cell_vec = sign(union_sig_cell_vec);
        union_sig_cells = find(union_sig_cell_vec == 1);
        
        %computing mean response trace for each simple odor stimulus presented
        resp_trace_mat = [];
        bad_cells_all = [];
        for stim_type_n = 1:size(sig_cell_mat_key, 1)
            curr_stim_vec = sig_cell_mat_key(stim_type_n, :);
            
            olf1_od_n = curr_stim_vec(1, 1);
            olf1_dur = curr_stim_vec(1, 2);
            olf2_od_n = curr_stim_vec(1, 3);
            olf2_dur = curr_stim_vec(1, 4);
            
            if olf1_dur ~= 0 && olf2_dur ~= 0
                curr_stim_type = 2;     %handover stimuli
            elseif olf2_dur == 0
                curr_stim_type = 0;     %simple stimulus, delivered on olf1
            elseif olf1_dur == 0
                curr_stim_type = 1;     %simple stimulus, delivered on olf2                
            end
            
            curr_sig_cells = find(sig_cell_mat(:, stim_type_n) == 1);
            curr_trs = find(stim_mat_simple(:, od_col_ns(1)) == olf1_od_n & stim_mat_simple(:, dur_col_ns(1)) == olf1_dur &... 
                                stim_mat_simple(:, od_col_ns(2)) == olf2_od_n & stim_mat_simple(:, dur_col_ns(2)) == olf2_dur);
            
            stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);  %computing stimulus on and off frame numbers for olf1 and olf2
            
            curr_mean_traces = squeeze(mean(dff_data_mat_f(:, union_sig_cells, curr_trs), 3, 'omitnan'));
%             del = curr_mean_traces < 0;
%             curr_mean_traces(del) = 0;
            
%             curr_sds = std(curr_mean_traces((stim_frs(1) - round(10./frame_time)):(stim_frs(1) - 2), :), [], 1, 'omitnan');
%             curr_means = mean(curr_mean_traces((stim_frs(1) - round(5./frame_time)):(stim_frs(1) - 2), :), 1, 'omitnan');
            


            %identifying mis-identified junk cells
            stim_on = stim_frs{1};
            stim_on = stim_on(1) - 1;
            bad_cells = mean(abs(curr_mean_traces(1:stim_on, :)), 1, 'omitnan');
            bad_cells = find(bad_cells > 0.10);
            bad_cells_all = [bad_cells_all, bad_cells];
            
            %logging current fly's response traces
            resp_trace_mat = pad_n_concatenate(resp_trace_mat, curr_mean_traces, 3, nan);
            
            
        end
        %getting rid of misidentified junk cells
        resp_trace_mat(:, bad_cells_all, :) = [];
        bad_cells_all = [];
        
        %pooling sig cells across flies
        all_resp_traces = pad_n_concatenate(all_resp_traces, resp_trace_mat, 2, nan);
        disp(curr_dir)
        
     end
     
     od_name_vec = [];
     q_resp_mat = [];
     for stim_type_n = 1:size(sig_cell_mat_key, 1)
        curr_stim_vec = sig_cell_mat_key(stim_type_n, :);

        olf1_od_n = curr_stim_vec(1, 1);
        olf1_dur = curr_stim_vec(1, 2);
        olf2_od_n = curr_stim_vec(1, 3);
        olf2_dur = curr_stim_vec(1, 4);

        if olf1_dur ~= 0 && olf2_dur ~= 0
            curr_stim_type = 2;     %handover stimuli
        elseif olf2_dur == 0
            curr_stim_type = 0;     %simple stimulus, delivered on olf1
        elseif olf1_dur == 0
            curr_stim_type = 1;     %simple stimulus, delivered on olf2                
        end
        
        
        %piggy-backing on last dataset to determine stim_frs, assuming they
        %are the same across all datasets for the current stimulus type
        curr_trs = find(stim_mat_simple(:, od_col_ns(1)) == olf1_od_n & stim_mat_simple(:, dur_col_ns(1)) == olf1_dur &... 
                                stim_mat_simple(:, od_col_ns(2)) == olf2_od_n & stim_mat_simple(:, dur_col_ns(2)) == olf2_dur);
            
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);  %computing stimulus on and off frame numbers for olf1 and olf2
        
        %constructing odor stimulus name string
        olf1_od_ind = find(odor_list_olf1 == olf1_od_n);
        olf1_od_name = od_names{olf1_od_ind};
        olf2_od_ind = find(odor_list_olf2 == olf2_od_n);
        olf2_od_name = od_names{olf2_od_ind};

        if curr_stim_type == 0  %simple stimulus on olf1
            curr_stim_name = [olf1_od_name, ' - single pulse'];
            stim_frs = stim_frs{1};
            curr_color = color_vec(olf1_od_ind, :);
        elseif curr_stim_type == 1 %simple stimulus on olf2
            curr_stim_name = [olf2_od_name, ' - single pulse'];
            stim_frs = stim_frs{2};
            curr_color = color_vec(olf2_od_ind, :);
        elseif curr_stim_type == 2 %odor transition stimulus
            curr_stim_name = [olf1_od_name, ' - ', olf2_od_name, ' transition'];
            stim_frs = [stim_frs{1}; stim_frs{2}];
            curr_color = [color_vec(olf1_od_ind, :); color_vec(olf2_od_ind, :)];
        else
        end
        
        curr_mean_traces = squeeze(all_resp_traces(:, :, stim_type_n));
        
        %sorting cells by response size to first stim_type
        if stim_type_n == 1
            wt_vec = mean(curr_mean_traces(stim_frs(1):stim_frs(2), :), 1, 'omitnan');     %response sizes on stim_type 1
        else
        end
        %sorting by wt_vec
        curr_mean_traces = [wt_vec; curr_mean_traces]';
        curr_mean_traces = sortrows(curr_mean_traces)';
        curr_mean_traces(1, :) = [];
        
        %computing mean resp over odor period (pulse2 for transition pulses)
        if curr_stim_type == 1
            q_resp_vec = max(curr_mean_traces(stim_frs(1, 1):(stim_frs(1, 2) + round(2./frame_time)), :));      %simple stimuli
            od_name_vec = [od_name_vec; {[olf2_od_name, ' simple']}];
        elseif curr_stim_type == 2
            q_resp_vec = max(curr_mean_traces(stim_frs(1, 1):(stim_frs(2, 2) + round(2./frame_time)), :));      %transition stimuli
            od_name_vec = [od_name_vec; {[olf2_od_name, ' transition']}];
        else
        end
        q_resp_mat = [q_resp_mat; q_resp_vec];
       
        
        fig_h = figure('Name', curr_stim_name);
        imagesc(curr_mean_traces', [0, 1.5]);
        colormap(greymap);
        set_xlabels_time(fig_h, frame_time, 10);
        ylabel('cell number');
        fig_wrapup(fig_h, []);
        add_stim_bar(fig_h, stim_frs, curr_color);
    end
    
     
    %plotting scatter plots and computing corrcoefs
    [corr_mat, p_mat] = corrcoef(q_resp_mat');
    
    pair_list = [1, 2;...   %1. Simple PA - Simple BA
                 1, 3;...   %2. Simple PA - Simple EL
                 2, 3;...   %3. Simple BA - Simple EL
                 5, 4;...   %4. Transition PA - Transition BA
                 1, 5;...   %4. Simple PA - Transition PA
                 2, 4;...   %4. Simple BA - Transition BA
                 ];
    
    for pair_n = 1:size(pair_list, 1)
        curr_pair = pair_list(pair_n, :);
        fig_h = figure('Name', [od_name_vec{curr_pair(1)}, ' - ', od_name_vec{curr_pair(2)}]);
        plot(q_resp_mat(curr_pair(1), :), q_resp_mat(curr_pair(2), :), 'O', 'MarkerSize', 4, 'MarkerEdgeColor', [0.65, 0.65, 0.65])
        xlabel([od_name_vec{curr_pair(1)}, ' responses']);
        ylabel([od_name_vec{curr_pair(2)}, ' responses']);
        axis([-0.2, 3, -0.2, 3]);
        ax_vals = axis;
        x_val = ax_vals(1) + ax_vals(2).*0.05;
        y_val = ax_vals(4).*0.95;
        r_val = num2str(round(corr_mat(curr_pair(1), curr_pair(2)), 3));
        p_val = num2str(round(p_mat(curr_pair(1), curr_pair(2)), 3));
        if p_val == '0'
            p_val = '< 0.001';
        else
        end
        text(x_val, y_val, ['r ', r_val, ', p ', p_val])
        fig_wrapup(fig_h, []);
    end
    
    
end

