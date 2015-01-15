function [fit_params, Yhat, se, ci] = linear_fit_stats(X, Y)
% function [fit_params, Yhat, se, ci] = linear_fit_stats(X, Y)
%
% This function returns the linear fit parameters to the data, along with
% a standard error and 95% condifence intervals.  X and Y must be vectors
% of the same size.

%let's eliminate NaNs right up front.  Easier than requiring the calling function
% to do it.
nnx = ~isnan(X); nny = ~isnan(Y); 
nn = nnx & nny; %use the vector and
X=X(nn); Y=Y(nn); 

n = length(X);
fit_params = polyfit(X,Y,1); %fit_params(1) is b, fit_params(2) is a in Y=a+bX
Yhat = polyval(fit_params, X);
residuals = Y-Yhat;
s2 = sum(residuals.^2)./(n-2);
s = sqrt(s2);
Xbar = mean(X);
se(2) = s*(sqrt((1/n) + (Xbar.^2 ./ sum((X-Xbar).^2)))); %se(a)
se(1) = s / sqrt(sum((X-Xbar).^2)); %se(b)
se = se';
fit_params = fit_params';
% 95% confidence intervals for fit parameters
ci = tinv(.975, n-2).*se;

if (1==0)
   figure; plot(X, Y, '.');
   figure; plot(X, residuals); 
end

stats.Yhat = Yhat;
stats.se = se;
stats.ci = ci;


