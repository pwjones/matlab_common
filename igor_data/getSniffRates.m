function [meanRates, allRates, respTrace, sniffs] = getSniffRates(exp, trialNum, frames, plotb)
% function sniffRates = getSniffRates(exp, trialNum, frames, plotb)
%

if isempty(frames)
    frames = 1:exp.vids(trialNum).nFrames;
end
fseg = listToIntervals(frames);
nseg = size(fseg,1);
freq = []; allRates = []; respTrace = []; sniffs = [];
for jj = 1:nseg
    frameInds = exp.camTrig(trialNum).frameInds(fseg(jj,1)):exp.camTrig(trialNum).frameInds(fseg(jj,2));
    freq = cat(1, freq, exp.resp(trialNum).sniffFreq(frameInds));
    mean_seg_freq(jj) = mean(exp.resp(trialNum).sniffFreq(frameInds));
    respTrace = cat(1, exp.resp(trialNum).value(frameInds));
    sniffs = cat(1, sniffs, exp.resp(trialNum).sniffVect(frameInds));
    %allRates = cat(1, allRates, freq(:));
end
meanRates = nanmean(freq);
allRates = freq;
if plotb % Then produce figures plotting the respiration frequencies over the tracked positions
    exp.vids(trialNum).plotPosition(frames, [], 0, 'k', '.');
    pos = exp.resp(trialNum).sniffPos;
    sniff_frames = exp.resp(trialNum).sniffFrames(exp.resp(trialNum).vidSniffs);
    [sniff_frames, sfi,~] = intersect(sniff_frames, frames);
    pfreq = exp.resp(trialNum).sniffFreq(exp.camTrig(trialNum).frameInds(sniff_frames)); 
    pos = pos(sfi, :);
    [cm cinds] = getIndexedColors('jet', pfreq);
    cm = brighten(cm, .7);
    for jj=1:length(pfreq)
        line('Parent',gca,'Xdata',pos(jj,1),'Ydata',pos(jj,2),'Marker','.','MarkerSize',8,'Color',cm(cinds(jj),:));
    end
    % to get the scale
    figure; colormap(cm);
    pcolor([pfreq, pfreq]);
    colorbar;
end
