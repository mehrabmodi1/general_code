function [J, grad] = costFunctionReg(theta, X, y, lambda)
%COSTFUNCTIONREG Compute cost and gradient for logistic regression with regularization
%   J = COSTFUNCTIONREG(theta, X, y, lambda) computes the cost of using
%   theta as the parameter for regularized logistic regression and the
%   gradient of the cost w.r.t. to the parameters. 

% Initialize some useful values
m = length(y); % number of training examples

% You need to return the following variables correctly 
J = 0;
grad = zeros(size(theta));

% ====================== YOUR CODE HERE ======================
% Instructions: Compute the cost of a particular choice of theta.
%               You should set J to the cost.
%               Compute the partial derivatives and set grad to the partial
%               derivatives of the cost w.r.t. each parameter in theta


y_pred1 = X*theta;
h_theta = sigmoid(y_pred1);

J1 = 1./m.*(sum(-y.*log(h_theta) - (1 - y).*log(1 - h_theta)));
Jlambda = lambda./(2.*m).*sum(theta(2:end).^2);     %computing regularisation component of cost excluding theta(1)
J = J1 + Jlambda;
grad1 = 1./m.*((h_theta - y)'*X);
grad1 = grad1(1);   %computing gradient for theta(1) without regularisation
grad = 1./m.*((h_theta - y)'*X) + (lambda./m).*theta';
grad(1) = grad1;    %subbing in un-regularised gradient for theta(1)




% =============================================================

end
