function [mu, ci] = eCI(data, alpha)
% function [mu, ci] = eCI(data)
% 
% Empirical Confidence intervals
% Data should be given in columns
% Alpha gives the significance level of the interval (eg .95)


mu = nanmean(data);
ci = zeros(2,size(data,2));
for ii = 1:size(data,2)
    [f,xi] = ksdensity(data(:,ii) ,'function', 'cdf');
    lb = find(f <= (1-alpha), 1, 'last');
    ci(1,ii) = xi(lb);
    ub = find(f >= alpha, 1, 'first');
    ci(2,ii) = xi(ub);
end
