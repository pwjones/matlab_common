function filt_trace = gaussianFilter(in_trace, stdg, varargin)
% function filt_trace = gaussianFilter(in_trace, stdg)
%
% Function to guassian filter a single vector with a kernel of
% the standard deviation given by stdg (in number of samples). 
% Varargin gives the choice of function to filter with: conv2 or filtfilt.
% Default is filtfilt.

if stdg ~= 0
    if ~isempty(varargin) && strcmp(varargin{1}, 'conv');
        convb = 1;
    else
        convb = 0;
    end
    
    filtx = -3*stdg:1:3*stdg;
    filty = normpdf(filtx,0,stdg);
    filty = filty/sum(filty);

    if convb %there is a better way to do this, but since they have different args, this is easy
        filt_trace = conv2(in_trace(:), filty(:), 'same');
    else
        filt_trace = filtfilt(filty, 1, in_trace);
    end
    
    x = 1:length(in_trace);
    %figure; plot(x, in_trace, x, filt_trace);
else
    filt_trace = in_trace;
end