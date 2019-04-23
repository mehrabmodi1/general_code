function [param_mat] = setup_param_mat_modular(params, stimGridFields)
%This function is called by setUpStimuli_modular to create a stimulus
%parameter matrix that includes combinations of various stim params for a
%given olfactometer.


for pp=1:length(params) %this loop is for each parameter, not trials
    
    %looping through each stim param of interest to build a large cell array of
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


end