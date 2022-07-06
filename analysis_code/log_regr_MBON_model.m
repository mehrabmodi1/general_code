function [theta, h_theta, test_acc_simp, test_acc_transition, h_theta1] = log_regr_MBON_model(inputs_train, outputs_train, inputs_test, outputs_test, initial_wts_offset, lambda)
%This function is custom-written to train a logistic regressor on simple
%pulse KC response vectors with analog, simple pulse MBON responses as the
%output. The goal is to then decode transition KC responses with the
%trained model.

%re-arranging odor dimension as more trials in input and output arrays
X = [];
y = [];
X_test = [];
y_test = [];

for od_n = 1:size(inputs_train, 3)
    X = pad_n_concatenate(X, squeeze(inputs_train(:, :, od_n)), 2, nan);
    y = pad_n_concatenate(y, squeeze(outputs_train(:, od_n)), 1, nan);
end

for od_n = 1:size(inputs_test, 3)
    X_test = pad_n_concatenate(X_test, squeeze(inputs_test(:, :, od_n)), 2, nan);
end

X = X';
X_test = X_test';

%adding offset unit (all ones)
X = [zeros(size(X, 1), 1), X];
X(:, 1) = X(:, 1) + 1;
X_test = [zeros(size(X_test, 1), 1), X_test];
X_test(:, 1) = X_test(:, 1) + 1;

% %getting rid of nan rows in X
% nanrows = sum(isnan(X), 2);
% nanrows = find(nanrows > 0);
% X(nanrows, :) = [];
% y(nanrows) = [];

initial_theta = rand(size(X, 2), 1).*initial_wts_offset;  %initializing weights as random numbers

%setting up and running fit function
% Set Options
lambda = lambda ;        %user-specified regularization constant
options = optimoptions(@fminunc,'Algorithm','Quasi-Newton','GradObj', 'on', 'MaxIter', 1000);

%Optimize
try
    [theta, J, exit_flag] = fminunc(@(t)(costFunctionReg(t, X, y, lambda)), initial_theta, options);
catch
    theta = [];
    h_theta = [];
    
end


%fitting logistic regressor
% try
%     [B,dev,stats] = mnrfit(X, y);
% catch
%     keyboard
% end


%decoding single pulse test trial
y_pred = X_test*theta;
h_theta = sigmoid(y_pred);
test_acc_simp = mean(round(h_theta(1:3))' == outputs_test(1:3));    %fraction of left out test outputs correct
test_acc_transition = mean(round(h_theta(4:5))' == outputs_test(4:5));
theta1 = ones(size(theta, 1), size(theta, 2));
theta(1) = [];

%computing regressor output with all weights set to 1 (to measure pre-learning performance)
y_pred1 = X_test*theta1;
h_theta1 = sigmoid(y_pred1);



