function [intervals, desc] = listToIntervals(list)
% function intervals = listToIntervals(list)
%
% This function breaks up a possibly discontinuous list
% into a set of bounded intervals, and a flag saying if
% it was continuous or not.

% regularize list
list = sort(list);
list = unique(list);

diffs = diff(list);
breaks = find(diffs > 1);
intervals = zeros(length(breaks)+1, 2);
intervals(1,1) = list(1);
for i=1:(length(breaks))
    if ~isempty(breaks) %also will be only thing that executes in only loop
        intervals(i,2) = list(breaks(i));
        intervals(i+1, 1) = list(breaks(i)+1);
    end
end
intervals(length(breaks)+1, 2) = list(end);
if (size(intervals) > 1)
    desc = 'discontinuous';
else
    desc = 'continuous';
end
end