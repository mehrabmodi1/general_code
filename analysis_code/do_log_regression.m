function [decoder_resp_traces, pred_accuracy_mat] = do_log_regression(single_traces_mat, n_reps_vec, sig_cell_mat_key, train_stim_frs)
%This function is custom-written to train a logistic regressor for each
%stim type and 
frame_time = 0.099;

if size(train_stim_frs, 1) > 1
    error('expected stim_frs pre-selected for only 1 olfactometer.')
else
end

stim_win = train_stim_frs;
stim_win2 = train_stim_frs + diff(train_stim_frs);      %for transition trials, the second pulse is 5s later

 n_single_pulse_reps = sum(n_reps_vec(1:3));

%computing input vectors for training as mean resps during the odor window relative to stim_frs for every trial
input_mat_sing = squeeze(mean(single_traces_mat(stim_win(1):stim_win(2), :, 1:n_single_pulse_reps), 1, 'omitnan'));
input_mat_trans = squeeze(mean(single_traces_mat(stim_win2(1):stim_win2(2), :, (n_single_pulse_reps + 1):end), 1, 'omitnan'));
input_mat = [input_mat_sing, input_mat_trans];

%looping through each stim type and training a decoder for it
last_used_tr = 0;
theta_fit_mat = zeros(size(input_mat, 1), size(sig_cell_mat_key, 1));
prob_vecs_all = [];
acc_vecs_all = [];
for stim_type = 1:size(sig_cell_mat_key, 1)
    curr_trs = [(last_used_tr + 1):1:(last_used_tr + n_reps_vec(stim_type, 1))];     %set of traces in single_traces_mat that corresspond to current stimulus type
    last_used_tr = curr_trs(end);
    
    y_vec = zeros(size(single_traces_mat, 3), 1);        
    y_vec(curr_trs) = 1;                                            %vector of correct answers for current decoder
    
   
    %Fitting single and transition trials to trials only of their own type
    if stim_type <= 3
        curr_trs_all =  1:n_single_pulse_reps;    %list of all single pulse trials
    elseif stim_type > 3
        curr_trs_all = (n_single_pulse_reps + 1):length(y_vec);    %list of all single pulse trials
    else
    end
    
    %doing all possible leave-one-out fits and decoding all left-out trials to assess performance
    acc_vec = [];
    prob_vecs_curr = [];
    for l_tr_n = 1:length(curr_trs_all)
        
        curr_trs_train = curr_trs_all;
        curr_trs_train(l_tr_n) = [];
        
        
        %listing relevant variables and renaming for convenience
        X = input_mat(:, curr_trs_train)';  %pop-vectors of response sizes for all trials arranged as m x n, where m is the number of training examples and n the number of features (neurons)
        
        %getting rid of nan rows in X
        nanrows = sum(isnan(X), 2);
        nanrows = find(nanrows > 0);
        X(nanrows, :) = [];
        
        y = y_vec(curr_trs_train);                %list of correct answers for current stimulus type (size n x 1)
        y(nanrows) = [];
        
        %since only the last two sets of trials are usd for transition decoders, need to account for that in trial index
        if stim_type > 3
            l_tr_n_adj = l_tr_n + n_single_pulse_reps;
        else
            l_tr_n_adj = l_tr_n;
        end
        
        X_test = input_mat(:, l_tr_n_adj)';
        y_test = y_vec(l_tr_n_adj);
        X_test_time = squeeze(single_traces_mat(:, :, l_tr_n_adj));
        
        initial_theta = zeros(size(X, 2), 1);  %initializing weights as random numbers

        
        %setting up and running fit function
        % Set Options
        lambda = 1 ;        %arbitrarily chosen regularization constant
        options = optimoptions(@fminunc,'Algorithm','Quasi-Newton','GradObj', 'on', 'MaxIter', 1000);
        
       
        
        % Optimize
        try
            [theta, J, exit_flag] = fminunc(@(t)(costFunctionReg(t, X, y, lambda)), initial_theta, options);
        catch
            keyboard
        end
        
        %decoding test trial
        y_pred = X_test*theta;
        h_theta = round(sigmoid(y_pred));
        if h_theta == y_test
            acc_vec = [acc_vec; 1];
        else
            acc_vec = [acc_vec; 0];
        end
        
        %decoding in time
        prob_vec = zeros(size(X_test_time, 1), 1);
        for frame_n = 1:size(X_test_time, 1)
            curr_X = squeeze(X_test_time(frame_n, :));
            y_pred = curr_X*theta;
            h_theta = sigmoid(y_pred);
            prob_vec(frame_n, 1) = h_theta;
        end
        %logging time-decoded vecs for each left-out trial
        prob_vecs_curr = pad_n_concatenate(prob_vecs_curr, prob_vec, 2, nan);   %logging across left-out repeats of current stim_type
       
    end
    prob_vecs_all = pad_n_concatenate(prob_vecs_all, prob_vecs_curr, 3, nan);    %logging across stim_type
    acc_vecs_all = pad_n_concatenate(acc_vecs_all, acc_vec, 2, nan);
end
decoder_resp_traces = prob_vecs_all;
pred_accuracy_mat = acc_vecs_all;






