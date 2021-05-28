function [stim_mat, stim_mat_simple, column_heads, color_vec, good_tr_list, params_orig, PID_traces_matched, PID_traces_orig] = load_params_trains_modular(direc, tif_datenums, frame_time)
%syntax: [stim_mat, stim_mat_simple, column_heads, color_vec, bad_tr_list] = load_params_trains(direc, tif_datenums, match_tifs)
%This function compares the time stamps of the tiff files in a dataset with
%those saved for each trial in the params file and aligns the two sets of
%trial numbers based on these time stamps. It then extracts trial relevant
%stimulus data from the params file and saves it as a readable struct.

no_tifs = 0;
if isempty(tif_datenums) == 0
    datenum_check = 1;
    n_trials_t = size(tif_datenums, 2);
else
    disp('no tifs datenums, continuing to load params without matching trial times to tifs.');
    no_tifs = 1;
    datenum_check = 0;
    
    %reading in manually determined list of bad (z-drifted) trials
    try
        good_tr_list = load([direc, 'bad_trial_list.mat']);
        good_tr_list = good_tr_list.bad_tr_list;
    catch
        disp('no z-drift tr list found.')
        %keyboard
        good_tr_list = [];
    end
    

end

%identifying newest params file and loading it
prev_direc = pwd;
cd(direc)
dir_contents = dir('params*.mat' );
n_files = size(dir_contents, 1);
date_nums = zeros(n_files, 1) + nan;
for par_file_n = 1:n_files
    date_nums(par_file_n) = dir_contents(par_file_n).datenum;
end

[del, last_filen] = max(date_nums);
try
    par_filename = dir_contents(last_filen).name;
catch
    keyboard
end
params = load([direc, par_filename]);
params = params.params_mat;
params_orig = params;

if datenum_check == 1
    n_trials_p = size(params, 2);        %n_trials according to the param file
    if no_tifs == 1
        n_trials_t = n_trials_p;
    else
    end    
    
    %reading in PID traces for all trials actually acquired (trs_done > 0)
    for tr_n = 1:n_trials_p
        if params(tr_n).trs_done == 0
            tr_n = tr_n - 1;
            break
        else
        end
    end
    
    try
        [PID_traces, del, n_PIDdata_dims] = get_PID_traces(direc, 1:tr_n, frame_time, 1);
    catch
        %NOTE: check if frame_time input was given to
        %load_params_trains_modular
        keyboard
    end
   
    PID_traces_acq = PID_traces(:, 1:n_PIDdata_dims:end);    
    LED_traces_acq = PID_traces(:, 2:n_PIDdata_dims:end);    
    
    if n_PIDdata_dims > 2
        elec_traces_acq = PID_traces(:, 3:n_PIDdata_dims:end);    
    else
        elec_traces_acq = LED_traces_acq.*nan;
    end
    
    PID_traces_orig = pad_n_concatenate(PID_traces_acq, LED_traces_acq, 3, nan);
    PID_traces_orig = pad_n_concatenate(PID_traces_acq, elec_traces_acq, 3, nan);
    
    %padding with nans for un-acquired traces
    PID_traces = zeros(size(PID_traces_acq, 1), n_trials_p) + nan;
    PID_traces(:, 1:tr_n) = PID_traces_acq;
    
    LED_traces = zeros(size(LED_traces_acq, 1), n_trials_p) + nan;
    LED_traces(:, 1:tr_n) = LED_traces_acq;
    
    elec_traces = zeros(size(elec_traces_acq, 1), n_trials_p) + nan;
    elec_traces(:, 1:tr_n) = elec_traces_acq;
        
    %reading in manually determined list of bad (z-drifted) trials
    good_tr_list = load([direc, 'bad_trial_list.mat']);
    good_tr_list = good_tr_list.bad_tr_list;
        
    n_trials = min([n_trials_t, n_trials_p]);     %number of matchable trials
    match_mat = zeros(n_trials_p, n_trials_t) + nan;

    for trial_n_p = 1:n_trials_p
        curr_datenum_p = params(trial_n_p).trs_done;
        curr_latency_p = params(trial_n_p).stimLatency;
        for trial_n_t = 1:n_trials_t
            curr_datenum_t = tif_datenums(trial_n_t).tstamp;
            if isempty(curr_datenum_t) == 1
                continue
            else
            end
            
            try
                curr_datenum_p = datetime(curr_datenum_p,'ConvertFrom','datenum');
            catch
                curr_datenum_p = datetime(curr_datenum_p);
            end
            %making sure there's no am-pm error
            if (curr_datenum_p - curr_datenum_t) < hours(11)
                match_mat(trial_n_p, trial_n_t) = seconds(curr_datenum_p - curr_datenum_t);      %calculating time elapsed from param time stamp to tif time stamp
            elseif (curr_datenum_p - curr_datenum_t) > hours(11)
                match_mat(trial_n_p, trial_n_t) = seconds(abs(curr_datenum_p - curr_datenum_t - hours(12)));      %calculating time elapsed from param time stamp to tif time stamp
            else
            end
            
        end   
    end
  
    
    pairing_tr_n_vec = find_pairing_tr(params);
    
    %finding matches for each param time_stamp
    saved_matches = [];
    for trial_n_p = 1:n_trials_p
        curr_vec = match_mat(trial_n_p, :);       %vec of time differences with all tiff timestamps for current param timestamp
        [delay, matched_t] = min(abs(curr_vec), [], 'omitnan');
       
        %checking if currently matched tif_n has already been matched to
        %another par_n
        if trial_n_p > 1
            dup_i = find(saved_matches(:, 2) == matched_t);
            if isempty(dup_i) == 0                  %duplicate match found
                dup_delay = saved_matches(dup_i, 1);
                if dup_delay > delay                %duplicate already in saved_matches is falsely matched
                    saved_matches(dup_i, :) = [];
                    saved_matches = [saved_matches; [delay, matched_t, trial_n_p]];
                    %continue
                elseif dup_delay < delay            %duplicate already in saved_matches is the correct match
                %continue
                else
                end
            else
                
                saved_matches = [saved_matches; [delay, matched_t, trial_n_p]];
                
            end
        else
            saved_matches = [saved_matches; [delay, matched_t, trial_n_p]];
        end
        
        %getting rid of saved_matches where the delay between tif and trial timestamps is > 5s.
        long_delay_trs = find(saved_matches(:, 1) > 5);
        for long_delay_n = 1:length(long_delay_trs)
            curr_delayed_tr = long_delay_trs(long_delay_n);
            if isnan(saved_matches(curr_delayed_tr, 2)) == 0
                saved_matches(curr_delayed_tr, :) = [];
            else
            end
        end        
        
        %accounting for case when no tif was acquired on the pairing trial
        %within the matching loop
        if isempty(intersect(trial_n_p, pairing_tr_n_vec)) == 0
            if saved_matches(size(saved_matches, 1), 3) ~= trial_n_p
                %case where current trial is led pairing trial, but no tiff was acquired
                saved_matches = [saved_matches; [delay, nan, trial_n_p]];
                
            else
            end
           
        else
        end
        
    end
    
    %syncing up match tr list with good_tr_list
    good_tr_list_bool = zeros(n_trials_t, 1);
    good_tr_list_bool(good_tr_list) = 1;
    matched_tr_list = saved_matches(:, 2);
    
    %accounting for cases where a pairing trial had no tif associated with
    %it by inserting an extra trial for it, and offsetting subsequent
    %trials by 1 position in the vector
    good_tr_list_booli = [];
    for matched_tr_n = 1:length(matched_tr_list)
        if isnan(matched_tr_list(matched_tr_n)) == 0
            good_tr_list_booli = [good_tr_list_booli; good_tr_list_bool(matched_tr_list(matched_tr_n))];
        elseif isnan(matched_tr_list(matched_tr_n)) == 1
            good_tr_list_booli = [good_tr_list_booli; 1];
        else
        end
    end
    good_tr_list = find(good_tr_list_booli == 1);
    
    tif_num = saved_matches(:, 2);
    par_num = saved_matches(:, 3);    
    n_matched_trials = length(par_num);

%     figure(1)
%     imagesc(match_mat)
%     xlabel('tif trial n');
%     ylabel('param trial n');
%     
%     figure(2)
%     plot(saved_matches, '.');
%     xlabel('param trial n')
%     ylabel('tif trial n');

   
else
    %allowing for cases where there are no tif_times associated with the loaded params file
    if no_tifs == 1
        n_trials_t = size(params, 2);
        
    else
    end    
    
    n_matched_trials = n_trials_t;
    par_num = 1:1:n_trials_t;
    tif_num = 1:1:n_trials_t;
    
    tr_n = n_trials_t;
    n_trials_p = n_trials_t;
    try
        [PID_traces, del, n_PIDdata_dims] = get_PID_traces(direc, 1:tr_n, frame_time, 1);
    catch
        %NOTE: check if frame_time input was given to
        %load_params_trains_modular
        keyboard
    end
   
    PID_traces_acq = PID_traces(:, 1:n_PIDdata_dims:end);    
    LED_traces_acq = PID_traces(:, 2:n_PIDdata_dims:end);    
    
    if n_PIDdata_dims > 2
        elec_traces_acq = PID_traces(:, 3:n_PIDdata_dims:end);    
    else
        elec_traces_acq = LED_traces_acq.*nan;
    end
    
    PID_traces_orig = pad_n_concatenate(PID_traces_acq, LED_traces_acq, 3, nan);
    PID_traces_orig = pad_n_concatenate(PID_traces_acq, elec_traces_acq, 3, nan);
    
    %padding with nans for un-acquired traces
    PID_traces = zeros(size(PID_traces_acq, 1), n_trials_p) + nan;
    PID_traces(:, 1:tr_n) = PID_traces_acq;
    
    LED_traces = zeros(size(LED_traces_acq, 1), n_trials_p) + nan;
    LED_traces(:, 1:tr_n) = LED_traces_acq;
    
    elec_traces = zeros(size(elec_traces_acq, 1), n_trials_p) + nan;
    elec_traces(:, 1:tr_n) = elec_traces_acq;
    
    
    
end

column_heads = [{'odor_n'}, {'duration'}, {'pulse_type'}, {'n_odor_pulses'}, {'inter_pulse_interval'}, {'rand_trains'}, {'mean_rand_pulse_dur'}, {'pulse_train_n'}, {'odour_olf2'}, {'duration_olf2'}, {'rel_stim_latency_olf2'}, {'pulse_type_olf2'}, {'n_od_pulses_olf2'}, {'inter_pulse_interval_olf2'}, {'rand_trains_olf2'}, {'mean_rand_pulse_dur_olf2'}, {'pulse_train_n_olf2'}, {'led_on'}, {'elec_on'}, {'isi'}, {'stim_latency'}, {'post_od_scan_dur'}, {'first_dilution'}, {'second_dilution'}, {'matched_tif_n'}];

%loop to read param values from params file into stim_mat
stim_mat_simple = zeros(n_matched_trials, 25) + nan;
for trial_n = 1:n_matched_trials
    curr_tr_p = par_num(trial_n);
    curr_tr_t = tif_num(trial_n);
    
    stim_mat(trial_n) = params(curr_tr_p);
    
    PID_traces_matched(:, trial_n) = PID_traces(:, curr_tr_p);
    LED_traces_matched(:, trial_n) = LED_traces(:, curr_tr_p);
   
    stim_mat_simple(trial_n, :) = [params(curr_tr_p).odours, params(curr_tr_p).duration, params(curr_tr_p).pulse_type, params(curr_tr_p).n_od_pulses, params(curr_tr_p).inter_pulse_interval, params(curr_tr_p).rand_trains, params(curr_tr_p).mean_rand_pulse_dur, params(curr_tr_p).pulse_train_n, ...
    params(curr_tr_p).odours_olf2, params(curr_tr_p).duration_olf2, params(curr_tr_p).rel_stimLatency_olf2, params(curr_tr_p).pulse_type_olf2, params(curr_tr_p).n_od_pulses_olf2, params(curr_tr_p).inter_pulse_interval_olf2, ...
    params(curr_tr_p).rand_trains_olf2, params(curr_tr_p).mean_rand_pulse_dur_olf2,  params(curr_tr_p).pulse_train_n_olf2, params(curr_tr_p).led_on, params(curr_tr_p).elec_on, params(curr_tr_p).isi, params(curr_tr_p).stimLatency, params(curr_tr_p).post_od_scan_dur, params(curr_tr_p).firstDilution, params(curr_tr_p).secondDilution, curr_tr_t];
end


PID_traces_matched = pad_n_concatenate(PID_traces_matched, LED_traces_matched, 3, nan);

%adding on match tiff_n
for trial_n = 1:n_matched_trials
    curr_tr_t = tif_num(trial_n);
    stim_mat(trial_n).matched_tif_n = curr_tr_t;
end

save([direc, '\stim_mat.mat'], 'stim_mat')

%setting up the standard color vector to use for plotting depending on how
%many odors there are in this dataset
n_odors = length(unique(stim_mat_simple(:, 1))); 
color_vec = cbrewer('qual', 'Dark2' ,n_odors, 'cubic');

end