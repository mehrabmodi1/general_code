function [params, params_spec]=setUpStimuli_modular(params)
% function params=setUpStimuli(params)
% This function sets up a stimulus parameter matrix where each trial's
% stimulus is defined as a train of pulses. Single pulse trials are trains of
% length 1. Each pulse train is defined by an external function, allowing
% flexible pulse-train definitions.
% Rob Campbell - October 2009, re-written by Mehrab Modi - 2016, 2019.


    %stimulus parameters (more can be added if they are numeric fields
    %in params)
    pulse_type = params.pulse_type;
    stimGridFields={'duration','stimLatency','odours','isi','firstDilution','secondDilution', 'n_od_pulses', 'inter_pulse_interval', 'post_od_scan_dur', 'mean_rand_pulse_dur', 'rand_trains', 'rand_train'};
    %olf2 params are not included here because that would inflexibly apply all combinations of all params across the two olfactometers. Now, each
    %combination of olf1-olf2 params needs to be separately set up through repaeat calls of setUpStimuli_modular.

    
    %making sure that the train numbers are listed individually in param_mat below.
    del = max([params.n_rand_trains, params.n_rand_trains_olf2]);
    n_rand_trains = max(del);
    if params.rand_trains == 0 && params.rand_trains_olf2 == 0
        n_rand_trains = 1;
    else
    end
    
    if params.rand_trains == 1
        params.rand_train = 1:n_rand_trains; 
    else
        params.rand_train = ones(n_rand_trains);
    end
    
    if params.rand_trains_olf2 == 1
        params.rand_train_olf2 = 1:n_rand_trains;
    else
        params.rand_train_olf2 = ones(n_rand_trains);
    end
    
    
    
    param_mat = setup_param_mat_modular(params, stimGridFields);
    %rand_train is a boolean that specifies whether to use a train of random duration odor pulses that add up to duration. 

    %n_rand_trains specifies how many different random trains are to be generated - these are treated as
    %separate stimuli for each odor and repeated reps number of times.

    %creating a copy of the rand_train column (for olf1) to keep track of rand train
    %number
    param_mat = [param_mat, param_mat(:, 12)];
    
    
    %checking if olf2 is being used in current set of trials
    if params.odours_olf2 == 0
        no_olf2 = 1;
    elseif params.duration_olf2 == 0
        no_olf2 = 1;
    elseif isempty(params.odours_olf2) == 1
        no_olf2 = 1;
    elseif isempty(params.duration_olf2) == 1
        no_olf2 = 1;
    else
        no_olf2 = 0;
    end
        
    if no_olf2 == 1 
        
        params.odours_olf2 = zeros(1, length(params.odours));
        params.duration_olf2 = zeros(1, length(params.duration));
        params.n_od_pulses_olf2 = zeros(1, length(params.n_od_pulses));
        params.inter_pulse_interval_olf2 = zeros(1, length(params.inter_pulse_interval));
        params.rand_trains_olf2 = zeros(1, length(params.rand_trains));
        params.mean_rand_pulse_dur_olf2 = zeros(1, length(params.mean_rand_pulse_dur));
    else
    end
   
    
    %creating a new param_mat for olf2 and appending as extra columns to the main param_mat.
    %checking for inconsistencies between param list lengths for olf1 and olf2
    if length(params.odours) ~= length(params.odours_olf2)
        error('length of params.odours not equal to length of params.odours_olf2.')
    elseif length(params.duration) ~= length(params.duration_olf2)
        error('length of params.duration not equal to length of params.duration_olf2.')
    elseif length(params.n_od_pulses) ~= length(params.n_od_pulses_olf2)
        error('length of params.n_od_pulses not equal to length of params.n_od_pulses_olf2.')
    elseif length(params.inter_pulse_interval) ~= length(params.inter_pulse_interval_olf2)
        error('length of params.inter_pulse_interval not equal to length of params.inter_pulse_interval_olf2.')
    elseif length(params.rand_trains) ~= length(params.rand_trains_olf2)
        error('params.n_rand_trains not equal to params.n_rand_trains_olf2.') 
     elseif length(params.mean_rand_pulse_dur) ~= length(params.mean_rand_pulse_dur_olf2)
        error('length of params.mean_rand_pulse_dur not equal to length of params.mean_rand_pulse_dur_olf2.') 
    else
    end

    stimGridFields_olf2 = {'duration_olf2','rel_stimLatency_olf2','odours_olf2', 'n_od_pulses_olf2', 'inter_pulse_interval_olf2', 'mean_rand_pulse_dur_olf2', 'rand_trains_olf2', 'rand_train_olf2'};
    param_mat_olf2 = setup_param_mat_modular(params, stimGridFields_olf2);
    param_mat = [param_mat, param_mat_olf2];
    param_mat = [param_mat, param_mat(:, 21)];
    
    
    
    %setting up repeats with randomisation if specified
    unit_param_mat = param_mat;
    n_reps = params.reps;
    size_repeat = size(unit_param_mat, 1);
    param_mat = [];
    for rep_n = 1:n_reps
       
        if params.randomize == 1

            t_order = randperm(size_repeat)';
            sub_param_mat = [t_order, unit_param_mat];
            sub_param_mat = sortrows(sub_param_mat);
            sub_param_mat(:, 1) = [];

        else
            sub_param_mat = unit_param_mat;
        end
        param_mat = cat(1, param_mat, sub_param_mat);
    end
    n_trials = size(param_mat, 1);
    
   
    
    %generating lists of rows or trials with the same pulse train (list of on-off timings) before converting param_mat to a cell array
    n_durations = length(params.duration);
    n_rand_trains = params.n_rand_trains;
    mean_rand_pulse_durs = params.mean_rand_pulse_dur;
    n_mean_rand_pulse_durs = length(mean_rand_pulse_durs);
    row_lists = {};
    for duration_n = 1:n_durations
        curr_duration = params.duration(duration_n);
        for rand_train_n = 1:n_rand_trains
            for mean_rand_dur_n = 1:n_mean_rand_pulse_durs
                mean_rand_dur_ni = mean_rand_pulse_durs(mean_rand_dur_n);
                curr_row_list = find(param_mat(:, 1) == curr_duration & param_mat(:, 13) == rand_train_n & param_mat(:, 10) == mean_rand_dur_ni);
                row_lists = [row_lists; {curr_row_list}];       %these are lists of trials that should have a given pulse timing train - random or simple
            end
        end
    end
    
    %generating lists of rows or trials with the same pulse train (list of on-off timings) before converting param_mat to a cell array
    n_durations = length(params.duration_olf2);
    n_rand_trains = params.n_rand_trains_olf2;
    mean_rand_pulse_durs = params.mean_rand_pulse_dur_olf2;
    n_mean_rand_pulse_durs = length(mean_rand_pulse_durs);
    row_lists_olf2 = {};
    for duration_n = 1:n_durations
        curr_duration = params.duration_olf2(duration_n);
        for rand_train_n = 1:n_rand_trains
            for mean_rand_dur_n = 1:n_mean_rand_pulse_durs
                mean_rand_dur_ni = mean_rand_pulse_durs(mean_rand_dur_n);
                curr_row_list = find(param_mat(:, 14) == curr_duration & param_mat(:, 21) == rand_train_n & param_mat(:, 19) == mean_rand_dur_ni);
                row_lists_olf2 = [row_lists_olf2; {curr_row_list}];       %these are lists of trials that should have a given pulse timing train - random or simple
            end
        end
    end
    param_mat_orig = param_mat;
    param_mat = num2cell(param_mat);                        %converting to cell array in order to be able to handle stimulus trains
     
    
    
    
    %generating random pulse trains or simple trains and inserting them into the param_mat cell array at the appropriate train number for olf1
    if params.rand_trains == 1                  %Case1: random pulse trains
        %generating and inserting rand trains for olf1
        n_durations = length(params.duration);
        n_rand_trains = params.n_rand_trains;
        mean_rand_pulse_durs = params.mean_rand_pulse_dur;
        n_mean_rand_pulse_durs = length(mean_rand_pulse_durs);    
        r_train_counter = 1;
        for duration_n = 1:n_durations
            curr_duration = params.duration(duration_n);
            for rand_train_n = 1:n_rand_trains
                for mean_rand_pulse_dur_n = 1:n_mean_rand_pulse_durs
                   mean_rand_pulse_dur_ni = mean_rand_pulse_durs(mean_rand_pulse_dur_n);
                   rand_train = rand_train_generator(curr_duration, mean_rand_pulse_dur_ni, 0.5, 20);       %generating pulse trains with exp mean 8s, min dur 0.1s, max dur 20s 
                   curr_trs = row_lists{r_train_counter};           %list of rows in param_mat into which current train needs to be inserted
                   param_mat(curr_trs, 12) = {rand_train};          %inserting random train into the appropriate rows of param_mat
                   param_mat(curr_trs, [7, 8]) = {nan};             %forcing simple train parameters to nans because this is a rand train
                   r_train_counter = r_train_counter + 1;
                end
            end
        end
        
    elseif params.rand_trains == 0                %for case with simple trains, inserting trains of a single pulse of length equal to duration
      
        for duration_n = 1:n_durations
            curr_duration = params.duration(duration_n);
            curr_trs = row_lists{duration_n};    %list of rows in param_mat into which current train needs to be inserted
            
            for tr_n = 1:length(curr_trs)
                tr_ni = curr_trs(tr_n);
                simple_pulse_train = [];
                curr_n_pulses = param_mat{tr_ni, 7};
                curr_int_pul_i = param_mat{tr_n, 8};
                for pulse_n = 1:curr_n_pulses
                    if pulse_n == 1
                        curr_pulse = [0, curr_duration]; 
                    elseif pulse_n > 1
                         if pulse_type == 0
                            curr_pulse = [curr_duration, curr_duration];
                        elseif pulse_type == 1
                            curr_pulse = [curr_int_pul_i, curr_duration];
                        else
                        end
                    else
                    end
                    simple_pulse_train = [simple_pulse_train; curr_pulse];
                end
                param_mat(curr_trs, 12) = {simple_pulse_train};                     %inserting simple train into the appropriate rows of param_mat
                param_mat(curr_trs, 10) = {nan};             %forcing simple train parameters to nans because this is a rand train
            end
        end
    end
    
    
    
    
    %generating random pulse trains or simple trains and inserting them into the param_mat cell array at the appropriate train number for olf2
    if params.rand_trains_olf2 == 1
        %generating and inserting rand trains for olf2
        n_durations = length(params.duration_olf2);
        n_rand_trains = params.n_rand_trains_olf2;
        mean_rand_pulse_durs = params.mean_rand_pulse_dur_olf2;
        n_mean_rand_pulse_durs = length(mean_rand_pulse_durs);    
        r_train_counter = 1;
        for duration_n = 1:n_durations
            curr_duration = params.duration(duration_n);
            for rand_train_n = 1:n_rand_trains
                for mean_rand_pulse_dur_n = 1:n_mean_rand_pulse_durs
                   mean_rand_pulse_dur_ni = mean_rand_pulse_durs(mean_rand_pulse_dur_n);
                   rand_train = rand_train_generator(curr_duration, mean_rand_pulse_dur_ni, 0.5, 20);       %generating pulse trains with exp mean 8s, min dur 0.1s, max dur 20s 
                   curr_trs = row_lists_olf2{r_train_counter};           %list of rows in param_mat into which current train needs to be inserted
                   param_mat(curr_trs, 21) = {rand_train};          %inserting random train into the appropriate rows of param_mat
                   param_mat(curr_trs, [17, 18]) = {nan};             %forcing simple train parameters to nans because this is a rand train
                   r_train_counter = r_train_counter + 1;
                end
            end
        end
        
    elseif  params.rand_trains_olf2 == 0                %for case with simple trains, inserting trains of a single pulse of length equal to duration
        for duration_n = 1:n_durations
            curr_duration = params.duration_olf2(duration_n);
            curr_trs = row_lists_olf2{duration_n};    %list of rows in param_mat into which current train needs to be inserted
            
            for tr_n = 1:length(curr_trs)
                tr_ni = curr_trs(tr_n);
                simple_pulse_train = [];
                curr_n_pulses = param_mat{tr_ni, 17};
                curr_int_pul_i = param_mat{tr_n, 18};
                for pulse_n = 1:curr_n_pulses
                    if pulse_n == 1
                        curr_pulse = [0, curr_duration]; 
                    elseif pulse_n > 1
                         if pulse_type == 0
                            curr_pulse = [curr_duration, curr_duration];
                        elseif pulse_type == 1
                            curr_pulse = [curr_int_pul_i, curr_duration];
                        else
                        end
                    else
                    end
                    simple_pulse_train = [simple_pulse_train; curr_pulse];
                end
                param_mat(curr_trs, 21) = {simple_pulse_train};                     %inserting simple train into the appropriate rows of param_mat
                param_mat(curr_trs, 19) = {nan};             %forcing simple train parameters to nans because this is a rand train
            end
        end
    end
    
    params_spec = params;
    clear params
    
    
    %forcing olf2 params to nan if no olf2 pulses are being delivered
    if no_olf2 == 1
        param_mat(:, 14:22) = {nan};        
    else
    end
    
    
    %filling params from param_mat matrix into params structure
    %olf1 params
    stimGridFields = [stimGridFields, 'rand_train_n'];          %adding this here to add an extra field to the param structure to keep track of the rand train number
    stimGridFields{1, 12} = 'pulse_train';
    stimGridFields{1, 13} = 'pulse_train_n';
    for trial_n = 1:size(param_mat, 1)
        for s_param_n = 1:length(stimGridFields)
            s_param = stimGridFields{s_param_n};
            params(trial_n).(s_param) = param_mat{trial_n, s_param_n};
            
            %forcing inter_pulse_interval = duration for 0.5 duty cycle pulse-type just for record keeping consistency
            if pulse_type == 0 && s_param_n == 7
                params(trial_n).inter_pulse_interval = params(trial_n).duration;
            else
            end
                        
        end
        params(trial_n).odourNames = params_spec.odourNames;
    end
    
    %olf2 params
    stimGridFields = [stimGridFields_olf2, 'rand_train_n_olf2'];          %adding this here to add an extra field to the param structure to keep track of the rand train number
    stimGridFields{1, 8} = 'pulse_train_olf2';
    stimGridFields{1, 9} = 'pulse_train_n_olf2';
    for trial_n = 1:size(param_mat, 1)
        for s_param_n = 1:length(stimGridFields)
            s_param = stimGridFields{s_param_n};
            params(trial_n).(s_param) = param_mat{trial_n, (s_param_n + 13)};
            
            %forcing inter_pulse_interval = duration for 0.5 duty cycle pulse-type just for record keeping consistency
            if pulse_type == 0 && s_param_n == 7
                params(trial_n).inter_pulse_interval = params(trial_n).duration;
            else
            end
                        
        end
        params(trial_n).odourNames_olf2 = params_spec.odourNames_olf2;
    end
    
    
    
    %creating new fields in params to keep track of completed trials and
    %whether multi-pulses are 50% duty-cycle.
    for trial_n = 1:n_trials
        params(trial_n).trs_done = 0;
        params(trial_n).pulse_type = params_spec.pulse_type;
        
        if no_olf2 == 0
            params(trial_n).pulse_type_olf2 = params_spec.pulse_type_olf2;
        elseif no_olf2 == 1
            params(trial_n).pulse_type_olf2 = nan;
        else
        end
    end
    
    
    %checking for LED and elec odours and creating new fields in the output params matrix for LED and ELEC stim.
    led_ods = params_spec.led_odours;
    elec_ods = params_spec.elec_odours;
    for tr_n = 1:size(params, 2)
        %LED
        if isempty(intersect(params(tr_n).odours, led_ods)) == 0
            params(tr_n).led_on = 1;
            try
                params(tr_n).led_power = params_spec.LED_power;
            catch
                params(tr_n).led_power = 0;
            end
        elseif isempty(intersect(params(tr_n).odours, led_ods)) == 1
            params(tr_n).led_on = 0;
            params(tr_n).led_power = 0;
        else
        end


        %elec
        if isempty(intersect(params(tr_n).odours, elec_ods)) == 0
            params(tr_n).elec_on = 1;
        elseif isempty(intersect(params(tr_n).odours, elec_ods)) == 1
            params(tr_n).elec_on = 0;
        else
        end
        
        %adding other led/elec stimulus related parameters which will only
        %be used if elec_on == 1 or led_on == 1
        params(tr_n).rel_stim_init_delay = params_spec.rel_stim_init_delay;
        params(tr_n).stim_dur = params_spec.stim_dur;
        params(tr_n).stim_freq = params_spec.stim_freq;
        params(tr_n).st_duty_cyc = params_spec.st_duty_cyc;
        
    end
   
end


