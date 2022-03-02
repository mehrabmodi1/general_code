clear all
close all

dataset_list_paths = [ ...
                         {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\d5HT1b_similar_od_handovers.xls'};...
                         {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\c305a_similar_od_handovers.xls'};...
                         {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\c739_similar_od_handovers.xls'};...
                      ];

[del, odor_names1] = xlsread('C:\Data\Code\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
[del, odor_names2] = xlsread('C:\Data\Code\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList_olf2.xls', 1);
odor_names2{3} = 'Butyl acetate';
od_names = [{'PA'}, {'BA'}, {'EL'}];

save_path_base = 'C:\Data\Data\Analysed_data\Analysis_results\KC_transition_logregr\';
           
paired_color = [0,136,55]./256;
unpaired_color = [166,219,160]./256;
EL_color = [123,50,148]./256;
PA_color = [0,136,55]./256;
BA_color = [166,219,160]./256;
mean_color = [0.8, 0.4, 0.4];

suppress_QC_plots = 1;
pool_flies_traces = 0;

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
pause_PCAs = 0;
single_fly_n = 5;   %G KCs - 5, A'\B' KCs - 5, A/B KCs - 3
subtract_pulse1_resps = 0;      %1 - subtracts the appropriate mean responses to single pulses from the transition response traces
re_extract_KC_data = 0;


%extracting, re-formatting and writing KC response data to disk
if re_extract_KC_data == 1
    for list_n = 1:size(dataset_list_paths, 1)
        curr_dir_list_path = dataset_list_paths{list_n, 1};
        curr_dir_list_path = update_list_path(curr_dir_list_path);
        [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
        n_dirs = size(dir_list, 1);
        dataset_list_name = findstr(curr_dir_list_path, 'list_');
        dataset_list_name = curr_dir_list_path((dataset_list_name + 5):(end - 4));

        %parsing path to identify current KC type
        slashi = findstr(curr_dir_list_path, '\');
        slashi = slashi(length(slashi));
        uscorei = findstr(curr_dir_list_path, '_');
        uscorei = uscorei - slashi;
        uscorei = uscorei(uscorei > 0);
        KC_type = curr_dir_list_path((slashi + 1):(slashi + uscorei - 1));


        %loop to go through all experiment datasets listed in list file
        all_resp_traces = [];
        sig_cell_mat_all = [];
        saved_resp_sizes_all = [];
        n_cells_all = [];
        all_fly_accuracy_mat = [];
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
            [stim_mat, stim_mat_simple, column_heads, color_vec, good_tr_list, params_orig, PID_traces_matched, PID_traces_orig] = load_params_trains_modular(curr_dir, tif_times, frame_time);    %reading in trial stimulus parameters after matching time stamps to F traces
            color_vec = [PA_color; BA_color; EL_color];

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

            %sanity-checking KC data
            [del, all_bad_trs] = cell_data_quality_control(dff_data_mat_f, stim_mat, stim_mat_simple, column_heads, sig_cell_mat, sig_cell_mat_key, suppress_QC_plots, frame_time);

            union_sig_cell_vec = sum(sig_cell_mat, 2);
            union_sig_cell_vec = sign(union_sig_cell_vec);
            union_sig_cells = find(union_sig_cell_vec == 1);
            %sum(sig_cell_mat)./size(sig_cell_mat, 1);
            union_nonsig_cells = 1:size(sig_cell_mat, 1);
            union_nonsig_cells(union_sig_cells) = [];

            sig_cell_mat_all = [sig_cell_mat_all; sig_cell_mat(:, 1:3)];

            %computing mean response trace for each simple odor stimulus presented
            resp_trace_mat = [];
            bad_cells_all = [];
            saved_resp_sizes = [];
            single_traces_mat = [];
            n_reps_vec = [];
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
                    curr_stim_type_simp = 0;
                elseif olf1_dur == 0
                    curr_stim_type = 1;     %simple stimulus, delivered on olf2
                    curr_stim_type_simp = 1;
                end

                curr_sig_cells = union_sig_cells;
                curr_trs = find(stim_mat_simple(:, od_col_ns(1)) == olf1_od_n & stim_mat_simple(:, dur_col_ns(1)) == olf1_dur &... 
                                    stim_mat_simple(:, od_col_ns(2)) == olf2_od_n & stim_mat_simple(:, dur_col_ns(2)) == olf2_dur);

                %logging traces sorted by stim type to use for logistic regression decoding.
                single_traces_mat = pad_n_concatenate(single_traces_mat, dff_data_mat_f(:, :, curr_trs), 4, nan); 

                stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);  %computing stimulus on and off frame numbers for olf1 and olf2

                %computing response sizes and logging only for sig cells
                stim_frs_simp = stim_frs{2};
                stim_frs_simp(2) = stim_frs_simp(2) + round(2./frame_time);
                mean_resps(union_nonsig_cells) = 0;
                saved_resp_sizes(:, stim_type_n) = mean_resps';


            end

            %training logistic regression decoders
            %calling logistic regression main_function
    %         stim_frs_plt = stim_frs{1};
    %         try
    %             [decoder_resp_traces, pred_accuracy_mat] = do_log_regression(single_traces_mat, n_reps_vec, sig_cell_mat_key, stim_frs_plt);        
    %         catch
    %             keyboard
    %         end
    %         decoder_resp_traces = movmean(decoder_resp_traces, 5, 1);

            fly_data.resp_traces = single_traces_mat;
            fly_data.resp_sizes = saved_resp_sizes;
            fly_data.stim_mat_key = sig_cell_mat_key;
            fly_data.stim_frs = stim_frs;
            fly_data.stim_mat = stim_mat_simple;
            fly_data.KC_type = KC_type;
            

            save_dir = [save_path_base, KC_type, '\fly', num2str(dir_n)];
            mkdir(save_dir);
            save([save_dir, '\fly_data.mat'], 'fly_data');
            %keyboard
        end

       



    end
else
end


%Reading in KC and MBON response data



%reading in MBON data
MBON_path = [save_path_base, 'MBON\'];
MBON_simp_data = load([MBON_path, 'simple_stim_data.mat']);
MBON_simp_data = MBON_simp_data.all_MBON_data;
MBON_transition_data = load([MBON_path, 'transition_stim_data.mat']);
MBON_transition_data = MBON_transition_data.all_MBON_data;

%reading in KC data, doing fit and logging results
KC_types = [{'d5HT1b'}, {'c305a'}, {'c739'}];
for KC_type_n = 1:3
    KC_type = KC_types{KC_type_n};
    n_flies = length(dir([save_path_base, KC_type])) - 2;
    for fly_n = 1:n_flies
        KC_data = load([save_path_base, KC_type, '\fly', num2str(fly_n), '\fly_data.mat']);
        KC_data = KC_data.fly_data;
        KC_resp_data = KC_Data.resp_sizes;
        
        
        keyboard
        
        
    end
    
end




%constructing odor stimulus name string
olf1_od_ind = find(odor_list_olf1 == olf1_od_n);
olf1_od_name = od_names{olf1_od_ind};
olf2_od_ind = find(odor_list_olf2 == olf2_od_n);
olf2_od_name = od_names{olf2_od_ind};

if curr_stim_type == 0  %simple stimulus on olf1
    curr_stim_name = [olf1_od_name, ' - single pulse'];
    stim_frs = stim_frs{1};z                      
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
