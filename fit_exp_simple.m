function [x_vec, fit_params] = fit_exp_simple(y_vals, x_step)
%syntax: [fitted_params] = fit_exp_simple(y_vals)
%This function converts an exponential sample points into a linear sample
%in log-space. It then fits a line in log-space and returns the fit
%parameters which can then be used in linear space to model the original
%exponential. This is done because linear fitting works better.

y = y_vals;
y_vals_l = log(y_vals);         %converting exponential points into linear points in log-space
x_vec = 0;
x = 0;
for step_n = 1:(length(y_vals) - 1)
    x = x + x_step;
    x_vec = [x_vec; x];
end
if size(y_vals, 1) == 1
    tbl = table(x_vec, y_vals_l');
else
    tbl = table(x_vec, y_vals_l);
end

modelfun = @(b,x) b(1).*x + b(2);  %fitting a simple line
beta0 = [0.1, 1]; % Guess values to start with.  Just make your best guess.
try
    mdl = fitnlm(tbl, modelfun, beta0);
catch
    keyboard
end
fit_params = mdl.Coefficients{:, 'Estimate'};

%fit_params(2) = exp(fit_params(2));     %converting proportionality paramter for use in linear space