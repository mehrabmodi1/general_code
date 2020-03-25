clear all
close all

path_base = 'C:\Data\Data\Analysed_data\data_sharing\KC_long_trace\';

n_flies = 6;

for fly_n = 1:n_flies
    an_path = [path_base, 'fly', num2str(fly_n), '\'];
    
    %reading in KC response params fitted with Herve's program
    fit_params = readNPY([an_path, 'fit_params.npy']);
    fit_tr_list = load([an_path, 'tr_list.mat']);       %this is the list of trials with only long dur, single od pulse used for the KC param fits
    fit_tr_list = fit_tr_list.tr_list;
        
    %reading in measured KC response data and stimulus params
    curr_dir = load([an_path, 'path_orig.mat']);
    curr_dir = curr_dir.curr_dir;
    curr_dir = manage_base_paths(curr_dir, 3);
    curr_dir = raw_direc_with_1(curr_dir);
    curr_dir = [curr_dir, '\'];
    tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
    tif_times = tif_times.time_stamps;
    [stim_mat, stim_mat_simple, column_heads, color_vec, good_tr_list, params_orig] = load_params_trains_modular(curr_dir, tif_times);    %reading in trial stimulus parameters after matching time stamps to F traces%[stim_mat, stim_mat_simple, column_heads, color_vec, g_tr_list] = load_params_trains(curr_dir, tif_times);    %reading in trial stimulus parameters after matching time stamps to F traces

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


    %loading extracted raw fluorescence data matrices written by raw_dff_extractor
    raw_data_mat = load([curr_dir 'extracted_raw_data_mat.mat']);
    raw_data_mat = raw_data_mat.raw_data_mat;           %raw F traces extracted from ROIs
    raw_data_mat_orig = raw_data_mat;
    tif_n_col_n = find_stim_mat_simple_col('matched_tif_n', column_heads);

    %inserting dummy trials (all nans) into raw_data_mat for pairing trials for
    %which no corress .tifs were acquired
    raw_data_mat = match_up_rawmat_matchedtrs(raw_data_mat, stim_mat_simple, tif_n_col_n);

    bad_tr_list = 1:1:size(raw_data_mat, 3);
    bad_tr_list(good_tr_list) = [];
    raw_data_mat(:, :, bad_tr_list) = nan;
    n_cells = size(raw_data_mat, 2);
    
    %calculating dF/F traces from raw data
    filt_time = .5;            %in s, the time window for boxcar filter for generating filtered traces
    [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, curr_dir);


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
    
    %identifying trials used to fit KC params and comparing fits to mean traces
    fit_trs_olf1 = find(stim_mat_simple(fit_tr_list, dur_col_ns(1)) > 0.1);
    fit_trs_olf1 = fit_tr_list(fit_trs_olf1);
    fit_trs_olf1_ods = unique(stim_mat_simple(fit_trs_olf1, od_col_ns(1)));
    fit_trs_olf1_dur = unique(stim_mat_simple(fit_trs_olf1, dur_col_ns(1)));

    %plotting olf1 fits
    for od_n = 1:length(fit_trs_olf1_ods)
        od_ni = fit_trs_olf1_ods(od_n);
        curr_trs = find(stim_mat_simple(fit_trs_olf1, od_col_ns(1)) == od_ni);      %list of trials with current 
        param_od_ns = [od_n, nan];      %odor indices for olf1 and olf2; ie not odor numbers, but their positions in odor lists.     
        [olf1_model_resps, olf2_model_resps] = get_fit_traces(curr_trs, stim_mat, fit_params, param_od_ns, frame_time, dff_data_mat_f);
    end   
    
    
    fit_trs_olf2 = find(stim_mat_simple(fit_tr_list, dur_col_ns(2)) > 0.1);
    fit_trs_olf2 = fit_tr_list(fit_trs_olf2);
    fit_trs_olf2_ods = unique(stim_mat_simple(fit_trs_olf2, od_col_ns(2)));
    fit_trs_olf2_dur = unique(stim_mat_simple(fit_trs_olf2, dur_col_ns(2)));

    
    keyboard
    
end

function [olf1_model_resps, olf2_model_resps] = get_fit_traces(curr_trs, stim_mat, fit_params, param_od_ns, frame_time, dff_data_mat)
    
    %modelling responses to each pulse in ofl1 train 
    n_frames = size(dff_data_mat, 1);
    olf1_model_resps = zeros(n_frames, size(fit_params, 2)) + nan;
    curr_train_olf1 = stim_mat(curr_trs(1)).pulse_train;
    if size(curr_train_olf1, 1) == 1 && max(curr_train_olf1) == 0.1
        curr_train_olf1 = nan;
    else
    end
    if isnan(curr_train_olf1) == 0
        curr_latency = stim_mat(curr_trs(1)).stimLatency;
        curr_params = squeeze(fit_params(param_od_ns(1), :, :));        %fit params for current odor, for all cells
        curr_traces = dff_data_mat(:, :, curr_trs);
        curr_resp_traces = squeeze(mean(curr_traces, 3, 'omitnan'));
        for cell_n = 1:size(curr_params, 1)
            curr_resp_trace = curr_resp_traces(:, cell_n);
            curr_cell_params = curr_params(cell_n, :);
            olf1_model_resps(:, cell_n) = fit_KC_response(curr_train_olf1, curr_cell_params, frame_time, n_frames, curr_latency, curr_resp_trace);            
            
        end
    
    else
    end
    
    %modelling responses to each pulse in olf2 train
    curr_train_olf2 = stim_mat(curr_trs(1)).pulse_train_olf2;
    olf2_model_resps = zeros(n_frames, size(fit_params, 2)) + nan;
    curr_train_olf2 = stim_mat(curr_trs(1)).pulse_train_olf2;
    
    if isnan(curr_train_olf2) == 0
        curr_latency = stim_mat(curr_trs(1)).stimLatency;
        curr_params = squeeze(fit_params(param_od_ns(2), :, :));        %fit params for current odor, for all cells
        curr_traces = dff_data_mat(:, :, curr_trs);
        curr_resp_traces = squeeze(mean(curr_traces, 3, 'omitnan'));
        for cell_n = 1:size(curr_params, 1)
            curr_resp_trace = curr_resp_traces(:, cell_n);
            curr_cell_params = curr_params(cell_n, :);
            olf1_model_resps(:, cell_n) = fit_KC_response(curr_train_olf2, curr_cell_params, frame_time, n_frames, curr_latency, curr_resp_trace);            
            
        end
    
    else
    end

    
    
    
    function model_trace = fit_KC_response(curr_train, curr_cell_params, frame_time, n_frames, stim_latency, curr_resp_trace)
        a0 = curr_cell_params(1);
        a1 = curr_cell_params(2) .* -1;
        a2 = curr_cell_params(3) ;
        a3 = curr_cell_params(4) .* -1;
        t0 = curr_cell_params(5);
        t1 = curr_cell_params(6);
        t2 = curr_cell_params(7);
        t3 = curr_cell_params(8);
        
        n_pulses = size(curr_train, 1);
        model_trace = zeros(n_frames, 1);
        t = ((1:1:n_frames) * frame_time)';
        for pulse_n = 1:n_pulses
            curr_model_trace = zeros(n_frames, 1);
            curr_pulse = curr_train(pulse_n, :);
            curr_train_time = sum(sum(curr_train(1:(curr_pulse - 1), :)));
            on_time = stim_latency + curr_train_time + curr_pulse(1, 1);                       %time from beginning of train to onset of current pulse
            off_time = on_time + curr_pulse(1, 2);   %time from beginning of train to off of current pulse
            
            %boolean vectors that define the three epochs of the responses
            %to each pulse
            domain0 = t <= on_time;
            domain1 = (t > on_time) & (t <= off_time);
            domain2 = t > off_time;
            
            
            %constructing model response trace generated by current pulse
            curr_model_trace = curr_model_trace + domain0.*(a0);
            curr_model_trace = curr_model_trace + domain1.*(a1 .* (1.0 - exp(-(t - on_time) / t0))...
                                                            + a2 .* (1.0 - exp(-(t - on_time) / t1))...
                                                            + a0);
            %computing a scalar adjustment factor                                                        
            end1 = (a1 * (1.0 - exp(-on_time / t0))...
                    + a2 * (1.0 - exp(-on_time / t1))...
                    + a0);                       
            
            %adding the odor off part of the modelled response    
            curr_model_trace = curr_model_trace + domain2.*(a3 .* exp(-(t - off_time) / t2)...
                                                            + (end1 - a3) .* exp(-(t - off_time) / t3));                                            
            
            plot(curr_resp_trace);
            hold on
            plot(curr_model_trace', 'r')
            
            close figure 1
            
        end
        model_trace = model_trace + curr_model_trace;
       
    end
    
    
    
end