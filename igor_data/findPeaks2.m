function [peakTimes, dirs] = findPeaks2(signal, thresh, deriv_thresh, sign, peak_window, varargin)
% function [peakTimes, dirs] = findPeaks2(signal, thresh, deriv_thresh, sign, peak_window, varargin)
% 
% signal - the signal you want to detect peaks in
% thresh - the threshold value for peak detection
% deriv_thresh - the threshold put on the derivative that needs to be crossed (excludes small fluctations
%                around a flat line. 0 if you don't want it to be used.  
% sign - the direction, 1 for positive, -1 for negative, of peaks
% peak_window - the number of samples smaller than which you can't have multiple peaks.
%
% Variable arguments:
% filter len - to median filter the given signal within the peak detection function. To disable, give
%               length 0
% Refine peaks flag - 0: picks the lower derivative for peak 1: searches for the signal peak near the
%                   derivative zero crossing (1 works best, default)
 
debug = 0;
if nargin > 5
    filtlen = varargin{1};
else
    filtlen = 10; 
end
if nargin > 6 
    refine_peaks_b = varargin{2};
else
    refine_peaks_b = 1;
end


signal = signal(:);
if sign == -1
    signal = -signal;
end
eventVects = zeros(1, 2*peak_window+1);

thresdw = thresh/2;

n_pts = length(signal);
peakTimes = zeros(size(signal));
event_x = -peak_window:peak_window;

deriv = cat(1, 0, diff(signal));
if filtlen > 1
    deriv_filt = medfilt1(deriv, filtlen);
else
    deriv_filt = deriv;
end
signs = prod([deriv_filt(1:end-1), deriv_filt(2:end)]');
crossing_i = find(signs < 0);
ismax = signal(crossing_i+1) >= thresh;
maxi = crossing_i(ismax); % these are the local maxima

%deriv_thresh = .08; % used to weed out those that are going pretty much parallel to the trail
if deriv_thresh > 0
    keep = false(length(maxi),1);
    for ii = 1:length(maxi)
        dx = mean(abs(deriv_filt(maxi(ii):(maxi(ii)+1)))); 
        if dx >= deriv_thresh
            keep(ii)=1;
        end
    end
    maxi = maxi(keep);
end 
% this section finds the real peak within the window: gets the actual peak
% and not the peak derivative.

refined_max = NaN*zeros(length(maxi), 1); % new vector for the maxima
if ~isempty(maxi)
    n = 0;
    max_ind = 1; 
    while ~isempty(max_ind)
        n = n+1;
        i = maxi(max_ind);
        wind_min = max(1, i-peak_window);
        wind_max = min(n_pts, i+peak_window);
        if refine_peaks_b %picking the one with the maximum signal
            [~, li] = max(signal(wind_min:wind_max));
            refined_max(n) = li+wind_min;
        else %picking the one with the lowest derivative
            inds = maxi(maxi>= wind_min & maxi <= wind_max);
            [~, mini] = min(deriv_filt(inds)); 
            refined_max(n) = inds(mini);
        end
        max_ind = find(maxi > wind_max, 1, 'first');
    end
    refined_max = refined_max(1:n);
end

dirs = NaN*zeros(length(refined_max), 1); % Positive or negative peaks?
for ii=1:length(refined_max)
    wind_min = max(1, refined_max(ii)-1);
    slope = signal(refined_max(ii)) - signal(wind_min);
    if (slope > 0)
        dirs(ii)=1;
    else
        dirs(ii) = -1;
    end
end

if debug
    figure; hold on;
    plot(signal, 'k');
    plot(refined_max, signal(refined_max), 'ro');
end

peakTimes(refined_max) = 1;


