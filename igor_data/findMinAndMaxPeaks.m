function [peakVect, allPeakInds, allDirs] = findMinAndMaxPeaks(signal, thresh, peak_window, varargin)
% function peakVect = findMinAndMaxPeaks(signal, thresh, peak_window, varargin)

% findPeaks2 will find the peaks in a single direction.  Just do it both
% ways and respect the peak_window to noise peaks

[posPeaks, posDirs] = findPeaks2(signal, thresh, .05, 1, peak_window, varargin{:});
[negPeaks, negDirs] = findPeaks2(signal, thresh, .05, -1, peak_window, varargin{:});

posPeakInds = find(posPeaks);
negPeakInds = find(negPeaks);
allPeakInds = cat(1, posPeakInds(:), negPeakInds(:));
allDirs = cat(1, posDirs(:), -negDirs(:));
allPeaks = find(posPeaks + negPeaks);
[allPeakInds, ui, ~] = unique(allPeakInds);
allDirs = allDirs(ui);

peakVect = false(size(signal));
peakVect(allPeaks) = 1;

