function filt_trace = gaussianFilter(in_trace, stdg, varargin)
% function filt_trace = gaussianFilter(in_trace, stdg)
%
% Function to guassian filter a single vector with a kernel of
% the standard deviation given by stdg (in number of samples). 
% Varargin gives the choice of function to filter with: conv2 or filtfilt.
% Default is filtfilt.

if stdg ~= 0
    if nargin >= 3 && strcmp(varargin{1}, 'conv');
        convb = 1;
    else
        convb = 0;
    end
    if nargin >= 4 
        len_flag = varargin{2};
    else
        len_flag = 'same';
    end
    
    half_len = ceil(stdg*3);
    filtx = [-half_len:0, 1:half_len];
    filty = normpdf(filtx,0,stdg);
    filty = filty/sum(filty);

    if length(in_trace) >= 3*length(filtx)
        if convb %there is a better way to do this, but since they have different args, this is easy
            filt_trace = conv2(in_trace(:), filty(:), len_flag);
            if length(filt_trace) ~= length(in_trace)
                offset = floor(length(filtx)/2);
                pad = NaN*zeros(offset,1);
                filt_trace = cat(1, pad, filt_trace, pad);
            end
            if length(filt_trace) < length(in_trace)
                filt_trace = cat(1, filt_trace, NaN);
                %filt_trace = filt_trace(1:length(in_trace));
            end
        else
            filt_trace = filtfilt(filty, 1, in_trace);
        end
    else
        disp('Unable to filter - input must be 3x longer than filter');
        filt_trace = in_trace;
    end
    x = 1:length(in_trace);
    %figure; plot(x, in_trace, x, filt_trace);
else
    filt_trace = in_trace;
end