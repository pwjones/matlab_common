function [Coef, CoefCI, e, fit] = fitARmodel(data, p)
% function [Coef, e] = fitARmodel(data, p)
% 
% This is a simple, understandable function to fit an AR(p) model to a
% dataset DATA.  The form of the AR(p) model is given by
% X(i) = C + B(1)X(i-1) + B(2)X(i-2) + ..... + B(p)X(i-p)
% This function fits via least squares, C and B, returning them in Coef.
% e is a vector of the model residuals.
%

% pb is a boolean for plotting the residuals
pb = 0;

width = p;
 
y = data((width+2:end));

%# Construct matrix
X = zeros(length(y), width);
for ii = 1:width
    if ii == 1
        X(:,ii) = data((width+1):end-1);
    else
        X(:,ii) = conv(data((width-ii+1):end-2), (1/ii)*ones(ii,1), 'valid');
    end 
end
X = [ones(length(y),1) X];

%# Perform OLS
[Coef, CoefCI, e] = regress(y, X);
fit = X*Coef;
disp(['Explained variance: ' num2str(100 * (var(data) - var(e))./var(data)) '%']);

%# Perform a durbin watson test on the residuals
[DWpVal, DWStat] = dwtest(e, X);
if DWpVal < 0.05; 
    disp('WARNING: residuals from regression appear to be serially correlated. Estimated coefficients may not be consistent'); 
end

if pb
    ve = var(e);
    me = mean(e);
    [hy, bins] = hist(e,100);
    ny = normpdf(bins, me,ve);
    figure;
    bar(bins, hy./sum(hy));
    hold on;
    plot(bins, ny);
end