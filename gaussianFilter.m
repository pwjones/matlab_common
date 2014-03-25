function filt_trace = gaussianFilter(in_trace, stdg)
% function filt_trace = gaussianFilter(in_trace, stdg)
%
% Function to guassian filter a single vector with a kernel of
% the standard deviation given by stdg (in number of samples)

filtx = -3*stdg:1:3*stdg;
filty = normpdf(filtx,0,stdg);
filty = filty/sum(filty);

filt_trace = filtfilt(filty, 1, in_trace);