function [stim_mat, stim_mat_simple, column_heads, color_vec, good_tr_list] = load_params_trains(direc, tif_datenums, match_tifs)
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
    good_tr_list = load([direc, 'bad_trial_list.mat']);
    good_tr_list = good_tr_list.bad_tr_list;

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
par_filename = dir_contents(last_filen).name;
params = load([direc, par_filename]);


try
    params = params.params_mat;
catch
    keyboard
end



if datenum_check == 1
    n_trials_p = size(params, 2);        %n_trials according to the param file
    if no_tifs == 1
        n_trials_t = n_trials_p;
    else
    end    
    
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
            curr_datenum_p = datetime(curr_datenum_p,'ConvertFrom','datenum');
            
            %making sure there's no am-pm error
            if (curr_datenum_p - curr_datenum_t) < hours(11)
                match_mat(trial_n_p, trial_n_t) = seconds(curr_datenum_p - curr_datenum_t);      %calculating time elapsed from param time stamp to tif time stamp
            elseif (curr_datenum_p - curr_datenum_t) > hours(11)
                match_mat(trial_n_p, trial_n_t) = seconds(abs(curr_datenum_p - curr_datenum_t - hours(12)));      %calculating time elapsed from param time stamp to tif time stamp
            else
            end
        end   
    end
  
    %finding matches for each param time_stampg
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
                    continue
                elseif dup_delay < delay            %duplicate already in saved_matches is the correct match
                continue
                else
                end
            else
                
                saved_matches = [saved_matches; [delay, matched_t, trial_n_p]];
                
            end
        else
            saved_matches = [saved_matches; [delay, matched_t, trial_n_p]];
        end
    end
    
    %syncing up match tr list with good_tr_list
    good_tr_list_bool = zeros(n_trials_t, 1);
    good_tr_list_bool(good_tr_list) = 1;
    matched_tr_list = saved_matches(:, 2);
    good_tr_list_bool = good_tr_list_bool(matched_tr_list);
    good_tr_list = find(good_tr_list_bool == 1);
    
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
    %allowing for there to be no tifs associated with the loaded params file
    if no_tifs == 1
        n_trials_t = size(params, 2);
        
    else
    end    
    
    n_matched_trials = n_trials_t;
    par_num = 1:1:n_trials_t;
    tif_num = 1:1:n_trials_t;
end

column_heads = 'matched_tif_n, odor_n, duration, isi, n_odor_pulses, inter_pulse_interval, stim_latency, first_dilution, second_dilution, post_od_scan_dur, rand_train_n, multi_pulse_train, led_on, elec_on';


%loop to read param values from params file into stim_mat
stim_mat_simple = zeros(n_matched_trials, 14) + nan;
for trial_n = 1:n_matched_trials
    curr_tr_p = par_num(trial_n);
    curr_tr_t = tif_num(trial_n);
    stim_mat(trial_n).matched_tif_n = curr_tr_t;
    stim_mat(trial_n).odor_n = params(curr_tr_p).odours;
    stim_mat(trial_n).odor_duration = params(curr_tr_p).duration;
    stim_mat(trial_n).isi = params(curr_tr_p).isi;
    stim_mat(trial_n).n_odor_pulses = params(curr_tr_p).n_od_pulses;
    stim_mat(trial_n).inter_pulse_interval = params(curr_tr_p).inter_pulse_interval;
    stim_mat(trial_n).stim_latency = params(curr_tr_p).stimLatency;
    stim_mat(trial_n).first_dilution = params(curr_tr_p).firstDilution;
    stim_mat(trial_n).second_dilution = params(curr_tr_p).secondDilution;
    stim_mat(trial_n).post_od_scan_dur = params(curr_tr_p).post_od_scan_dur;
    stim_mat(trial_n).rand_trains = params(curr_tr_p).rand_train;
    stim_mat(trial_n).rand_train_n = params(curr_tr_p).rand_train_n;
    
    try
        stim_mat(trial_n).led_on = params(curr_tr_p).led_on;
        stim_mat(trial_n).elec_on = params(curr_tr_p).elec_on;
        stim_mat(trial_n).stim_init_delay_ms = params(curr_tr_p).stim_init_delay_ms;
        stim_mat(trial_n).stim_dur = params(curr_tr_p).stim_dur;
        stim_mat(trial_n).stim_freq = params(curr_tr_p).stim_freq;
        stim_mat(trial_n).st_duty_cyc = params(curr_tr_p).st_duty_cyc;
    catch
    end
    %checking if current trial was a rand train trial or a single pulse
    %trial
    if size(params(curr_tr_p).rand_train, 1) > 1
        train_on = 1;
    else
        train_on = 0;
    end
    
    try
       
        stim_mat_simple(trial_n, :) = [curr_tr_t, params(curr_tr_p).odours, params(curr_tr_p).duration, params(curr_tr_p).isi,...
            params(curr_tr_p).n_od_pulses, params(curr_tr_p).inter_pulse_interval, params(curr_tr_p).stimLatency,...
            params(curr_tr_p).firstDilution, params(curr_tr_p).secondDilution, params(curr_tr_p).post_od_scan_dur, params(curr_tr_p).rand_train_n, train_on, params(curr_tr_p).led_on, params(curr_tr_p).elec_on];
    catch
        if trial_n == 1
            stim_mat_simple(:, 13:14) = [];
        else
        end
        stim_mat_simple(trial_n, :) = [curr_tr_t, params(curr_tr_p).odours, params(curr_tr_p).duration, params(curr_tr_p).isi,...
            params(curr_tr_p).n_od_pulses, params(curr_tr_p).inter_pulse_interval, params(curr_tr_p).stimLatency,...
            params(curr_tr_p).firstDilution, params(curr_tr_p).secondDilution, params(curr_tr_p).post_od_scan_dur, params(curr_tr_p).rand_train_n, train_on];
    end

end
save([direc, '\stim_mat.mat'], 'stim_mat')

%setting up the standard color vector to use for plotting depending on how
%many odors there are in this dataset
n_odors = length(unique(stim_mat_simple(:, 2)));
color_vec = cbrewer('qual', 'Dark2' ,n_odors, 'cubic');



end