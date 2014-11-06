function peakTimes = findPeaks2(signal, thresh, sign, peak_window)
% function [peakTimes, eventVects] = findPeaks(signal, thresh, sign, peak_window)
% 
% signal - the signal you want to detect peaks in
% thresh - the threshold value for peak detection
% sign - the direction, 1 for positive, -1 for negative, of peaks
% peak_window - the number of samples smaller than which you can't have multiple peaks.
debug = 0;

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
deriv_filt = medfilt1(deriv, 10);
signs = prod([deriv_filt(1:end-1), deriv_filt(2:end)]');
crossing_i = find(signs < 0);
ismax = find(signal(crossing_i+1) >= thresh);
maxi = crossing_i(ismax)+1; % these are the local maxima

refined_max = NaN*zeros(length(maxi), 1); % new vector for the maxima
if ~isempty(maxi)
    n = 0;
    max_ind = 1; 
    while ~isempty(max_ind)
        n = n+1;
        i = maxi(max_ind);
        wind_min = max(1, i-peak_window);
        wind_max = min(n_pts, i+peak_window);
        [~, li] = max(signal(wind_min:wind_max));
        refined_max(n) = li+wind_min;
        
        max_ind = find(maxi > wind_max, 1, 'first');
    end
    refined_max = refined_max(1:n);
end

if debug
    figure; hold on;
    plot(signal, 'k');
    plot(refined_max, signal(refined_max), 'ro');
end

peakTimes(refined_max) = 1;


