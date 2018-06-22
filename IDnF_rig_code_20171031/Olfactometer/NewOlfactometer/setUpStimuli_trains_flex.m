function [params, params_spec]=setUpStimuli_trains_flex(params)
% function params=setUpStimuli(params)
% This function sets up a stimulus parameter matrix where each trial's
% stimulus is defined as a train of pulses. Single pulse trials are trains of
% length 1.

% Rob Campbell - October 2009, re-written by Mehrab Modi - 2016.


%stimulus parameters (more can be added if they are numeric fields
%in params)
pulse_type = params.pulse_type;
stimGridFields={'duration','stimLatency','odours','isi','firstDilution','secondDilution', 'n_od_pulses', 'inter_pulse_interval', 'post_od_scan_dur', 'rand_train'};

del = params.n_rand_trains;
n_rand_trains = max(del);
params.rand_train = 1:n_rand_trains; %making sure that the train numbers are listed individually in param_mat below.

%rand_train is a boolean that specifies whether to use a train of random
%duration odor pulses that add up to duration. 

%n_rand_trains specifies how many different random trains are to be generated - these are treated as
%separate stimuli for each odor and repeated reps number of times.

for pp=1:length(params) %this loop is not for trials
    
    %looping through each stim param of interest to build a large matrix of
    %all combinations of parameters of interest
    param_mat = [];                 %eventually, each row of this matrix will be a trial, and each column a parameter value
    for s_param_n = 1:length(stimGridFields)
        s_param = stimGridFields{s_param_n};
        s_param_vals = params.(s_param);
        
        if ischar(s_param_vals) == 1
            s_param_vals = str2num(s_param_vals);
        else
        end
        
        if size(s_param_vals, 1) < size(s_param_vals, 2)
            s_param_vals = s_param_vals';                   %making sure param vals are arranged in a column rather than a row
        else
        end
        
        n_vals = size(s_param_vals, 1);
        cur_mat_siz = size(param_mat, 1);
        
        if s_param_n == 1                    %initialising param mat
            param_mat = s_param_vals;
        
        else                                 %case when current parameter is not the first one in the list
            
            %loop to go through each value this parameter takes
            param_mat_orig = param_mat;
            param_mat = [];
            for par_val_n = 1:n_vals
                pad_vec = repmat(s_param_vals(par_val_n), cur_mat_siz, 1);   %col vector of current parameter value repeated so as to fill out one entire set of param_mat rows made so far
                param_mat = [param_mat; [param_mat_orig, pad_vec] ];
            end
            
        end
        
    end
    
    %creating a copy of the rand_train column to keep track of rand train
    %number
    param_mat = [param_mat, param_mat(:, 10)];
    
    %house-keeping steps
    n_reps = params.reps;
    param_mat = repmat(param_mat, n_reps, 1);      %replicating existing param_mat by n_reps
    n_trials = size(param_mat, 1);
    
    %randomising trial order if needed
    if params.randomize == 1
        t_order = randperm(n_trials)';
        param_mat = [t_order, param_mat];
        param_mat = sortrows(param_mat);
        param_mat(:, 1) = [];

    else
    end
    
    %generating random pulse trains and inserting them into param_mat as per rand train number, if needed
    n_durations = length(params.duration);
    %generating lists of rows before converting param_mat to a cell array
    row_lists = {};
    for duration_n = 1:n_durations
        curr_duration = params.duration(duration_n);
        for rand_train_n = 1:n_rand_trains
            curr_row_list = find(param_mat(:, 1) == curr_duration & param_mat(:, 10) == rand_train_n);
            row_lists = [row_lists; {curr_row_list}];
        end
    end
    
    param_mat = num2cell(param_mat);                %converting to cell array in order to be able to handle stimulus trains
    r_train_counter = 1;
    
    if params.rand_trains == 1
        for duration_n = 1:n_durations
            curr_duration = params.duration(duration_n);
            for rand_train_n = 1:n_rand_trains
               rand_train = rand_train_generator(curr_duration, 8, 0.1, 20);       %generating pulse trains with exp mean 5s, min dur 0.1s, max dur 10s 
               curr_trs = row_lists{r_train_counter};           %list of rows in param_mat into which current train needs to be inserted
               param_mat(curr_trs, 10) = {rand_train};                             %inserting random train into the appropriate rows of param_mat
               r_train_counter = r_train_counter + 1;
            end
        end
        
    elseif params.rand_trains == 0                %for case without random trains, inserting trains of a single pulse of length equal to duration
        
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
                param_mat(curr_trs, 10) = {simple_pulse_train};                     %inserting random train into the appropriate rows of param_mat
                
            end
        end
    end
    
    params_spec = params;
    clear params
    %filling params from param_mat matrix into params structure
    stimGridFields = [stimGridFields, 'rand_train_n'];          %adding this here to add an extra field to the param structure to keep track of the rand train number
    for trial_n = 1:size(param_mat, 1)
        for s_param_n = 1:length(stimGridFields)
            s_param = stimGridFields{s_param_n};
            params(trial_n).(s_param) = param_mat{trial_n, s_param_n};
            
            %forcing inter_pulse_interval = duration for 0.5 duty cycle pulse-type
            if pulse_type == 0 && s_param_n == 7
                params(trial_n).inter_pulse_interval = params(trial_n).duration;
            else
            end
                        
        end
        params(trial_n).odourNames = params_spec.odourNames;
    end
    
    %creating new fields in params to keep track of completed trials and
    %whether multi-pulses are 50% duty-cycle.
    for trial_n = 1:n_trials
        params(trial_n).trs_done = 0;
        params(trial_n).pulse_type = pulse_type;
    end
    
    
    %checking for LED and elec odours and creating new fields in the output params matrix for LED and ELEC stim.
    
    
    
end
keyboard
end

