function [peakVect, allPeakInds, allDirs] = findMinAndMaxPeaks(signal, thresh, peak_window, varargin)
% function peakVect = findMinAndMaxPeaks(signal, thresh, peak_window, varargin)

% findPeaks2 will find the peaks in a single direction.  Just do it both
% ways and respect the peak_window to noise peaks

[posPeaks, posDirs] = findPeaks2(signal, thresh, .08, 1, peak_window, varargin{:});
[negPeaks, negDirs] = findPeaks2(signal, thresh, .08, -1, peak_window, varargin{:});

posPeakInds = find(posPeaks);
negPeakInds = find(negPeaks);
allPeakInds = cat(1, posPeakInds(:), negPeakInds(:));
allDirs = cat(1, posDirs(:), -negDirs(:));
allPeaks = find(posPeaks + negPeaks);
[allPeakInds, ui, ~] = unique(allPeakInds);
allDirs = allDirs(ui);

peakVect = false(size(signal));
peakVect(allPeaks) = 1;

%% This is a bit of code that I thought I needed, but didn't end up working out.  
% peakVect = NaN*zeros(length(allPeaks), 1); % new vector for the maxima
% 
% %This bit of code just walks through the identified positive and negative peaks, picks the absolute value max
% %within the window around that peak, then moves onto the next identified peak that is outside the window.
% if ~isempty(peakVect)
%     n = 0;
%     max_ind = 1; 
%     while ~isempty(max_ind)
%         n = n+ 1;
%         i = allPeaks(max_ind);
%         wind_min = max(1, i-peak_window);
%         wind_max = min(n_pts, i+peak_window);
%         % excludes the identified peaks that are within the window
%         [~, li] = max(abs(signal(wind_min:wind_max)));
%         peakVect(n) = li+wind_min;
%         max_ind = find(allPeaks > wind_max, 1, 'first');
%     end
%     peakVect = peakVect(1:n);
% end

