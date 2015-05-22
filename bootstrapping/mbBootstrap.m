function [meanStat, ci, stat] = mbBootstrap(n, func, data, blockLen, blockOffset)
% function [meanStat, ci] = mbBootstrap(n, func, data, blockLen, blockOffset)
%
% Moving block bootstrap - generates a mean and CI around that mean for the
% moving blocks of timeseries input.
stat = [];
[y, q] = segments(data(:), blockLen, blockOffset);
tempdata = bootrsmp(n,y,q,numel(data)); % generates the matrix of resampled data
for ii=1:n 
    temp = func(tempdata(:,ii));
    stat = cat(1,stat, temp(:)');
end

meanStat = mean(stat);
ci = zeros(2,size(stat,2));
for ii = 1:size(stat,2)
    [f,xi] = ksdensity(stat(:,ii) ,'function', 'cdf');
    lb = find(f <= .05, 1, 'last');
    ci(1,ii) = xi(lb);
    ub = find(f >= .95, 1, 'first');
    ci(2,ii) = xi(ub);
end