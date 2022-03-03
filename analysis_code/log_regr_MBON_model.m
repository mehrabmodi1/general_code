function [theta, h_theta, test_acc, h_theta_transition, acc_transition] = log_regr_MBON_model(inputs_train, outputs_train, inputs_test, inputs_tr_test, paired_od_index)
%This function is custom-written to train a logistic regressor on simple
%pulse KC response vectors with analog, simple pulse MBON responses as the
%output. The goal is to then decode transition KC responses with the
%trained model.

%re-arranging odor dimension as more trials in input and output arrays
X = [];
y = [];
X_test = [];
for od_n = 1:3
    X = pad_n_concatenate(X, squeeze(inputs_train(:, :, od_n)), 2, nan);
    X_test = pad_n_concatenate(X_test, squeeze(inputs_test(:, :, od_n)), 2, nan);
    y = pad_n_concatenate(y, squeeze(outputs_train(1, od_n, :)), 1, nan);
end

%re-arranging transition trial KC responses such that paired odor responses are first and unpaired odor responses are second
X_test_transition = [];
if paired_od_index == 1
    for od_n = 1:2
        X_test_transition = pad_n_concatenate(X_test_transition, squeeze(inputs_tr_test(:, :, od_n)), 2, nan);
    end
elseif paired_od_index == 2
    for od_n = 2:-1:1
        X_test_transition = pad_n_concatenate(X_test_transition, squeeze(inputs_tr_test(:, :, od_n)), 2, nan);
    end
else
end

X = X';
y = y;
X_test = X_test';

%normalizing input and output vectors
X(X > 5) = 5;
X = X./5;
max_y = max(max(y, [], 'omitnan'), [], 'omitnan');
y = y./max_y;

%getting rid of nan rows in X
nanrows = sum(isnan(X), 2);
nanrows = find(nanrows > 0);
X(nanrows, :) = [];
y(nanrows) = [];


initial_theta = zeros(size(X, 2), 1);  %initializing weights as random numbers

%setting up and running fit function
% Set Options
lambda = 1 ;        %arbitrarily chosen regularization constant
options = optimoptions(@fminunc,'Algorithm','Quasi-Newton','GradObj', 'on', 'MaxIter', 1000);

%Optimize
try
    [theta, J, exit_flag] = fminunc(@(t)(costFunctionReg(t, X, y, lambda)), initial_theta, options);
catch
    theta = [];
    h_theta = [];
    keyboard
    return
end


%fitting logistic regressor
% try
%     [B,dev,stats] = mnrfit(X, y);
% catch
%     keyboard
% end



%decoding single pulse test trial
y_pred = X_test*theta;
h_theta = round(sigmoid(y_pred));
test_acc = sum(h_theta == [0; 0; 1])./3;    %fraction of left out test outputs correct


%decoding transition pulse test trial
y_pred_transition = X_test_transition'*theta;
h_theta_transition = round(sigmoid(y_pred_transition));
expected_h_theta_transition = zeros(16, 1);
expected_h_theta_transition(9:16, 1) = 1;
acc_transition = sum(h_theta_transition == expected_h_theta_transition)./16;




