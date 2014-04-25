function [yf, J] = logistic(p, xd)
%function yf = logistic(params, xd)
%
% Gives the results of the logistic function or logistic distribution CDF.  p are the parameters,
% with p(1) = height, p(2) = steepness, and p(3) = horizontal position.  
xd = xd(:);
yf = p(1)./(1+exp((-xd+p(3)).*p(2)));

if nargout > 1
    j1 = 1./(1+exp(-p(2).*(xd-p(3)))); %partial derivative, parameter 1
    j2 = p(1).*(xd-p(3))*exp(-p(2).*(xd-p(3)))./((1+exp((-xd+p(3)).*p(2))).^2);
    j3 = p(1).*p(2).*exp((-xd+p(3)).*p(2))./((1+exp((-xd+p(3)).*p(2))).^2);
    J = [j1, j2, j3]';
end