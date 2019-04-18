function params=setUpStimuli_EP(params)
% function params=setUpStimuli(params)
%
% Convert the parameters structure (defaultHBOparams) into a format
% for presenting to the fly. Everything will be presented from a
% BrainWare-like (www.tdt.com/software/for_neuro.htm) "stimulus
% grid", so we specify all stim parameters for each trial. This
% makes it easy to present whatever stimulus we want on a trial by
% trial basis.
%
%
% Rob Campbell - October 2009, re-written by Mehrab Modi - 2016.


%stimulus parameters (more can be added if they are numeric fields
%in params)
stimGridFields={'duration','stimLatency','odours','isi','firstDilution','secondDilution', 'n_od_pulses', 'inter_pulse_interval', 'post_od_scan_dur', 'rand_train', 'n_rand_trains'};


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


    %house-keeping steps
    n_reps = params.reps;
    param_mat = repmat(param_mat, n_reps, 1);      %replicating existing param_mat by n_reps
    param_mat = sortrows(param_mat);                                    %sorting
    n_trials = size(param_mat, 1);
    
    %randomising trial order if needed
    if params.randomize == 1
        t_order = randperm(n_trials)';
        param_mat = [t_order, param_mat];
        param_mat = sortrows(param_mat);
        param_mat(:, 1) = [];

    else
    end
    

    %filling params from param_mat matrix into params structure
    for s_param_n = 1:length(stimGridFields)
        s_param = stimGridFields{s_param_n};
        params.(s_param) = param_mat(:, s_param_n);
    end
    
    %creating new fields in params to keep track of completed trials
    params.trs_done = zeros(n_trials, 1);
   
    
    
end

end

