function [sniffRates, all_freq, all_pos, all_dist]  = getSniffRateOnTrail(exp, trials, trailNumber, threshDist, plotb)
% function plotSniffRateOnTrail(exp, trials, trailNumber, dist_thresh)
%

if isempty(trials) %allow for the option of leaving the trial numbers out
    trials = 1:length(exp.vids);
end

all_freq = []; all_pos = []; all_dist = [];
for ii = trials
    fseg = exp.vids(ii).getFollowingSegments([],trailNumber,threshDist);
    nseg = size(fseg,1);
    freq = []; frames = []; 
    for jj = 1:nseg
        frames = cat(2, frames, fseg(jj,1):fseg(jj,2)); %build a vector of frames from the intervals
        frameInds = exp.camTrig(ii).frameInds(fseg(jj,1)):exp.camTrig(ii).frameInds(fseg(jj,2));
        freq = cat(1, freq, exp.resp(ii).sniffFreq(frameInds));
        mean_seg_freq(jj) = mean(exp.resp(ii).sniffFreq(frameInds));
        %all_seg_freq = cat(1,all_seg_freq, exp.resp(ii).sniffFreq(frameInds));
    end
    sniffRates(ii) = nanmean(freq);
    
    
    if plotb % Then produce figures plotting the respiration frequencies over the tracked positions
        fr = frames;
        %fr = 1:exp.vids(ii).nFrames;
        exp.vids(ii).plotPosition(fr, [], 0, 'k', '.');
        pos = exp.resp(ii).sniffPos;
        sniff_frames = exp.resp(ii).sniffFrames(exp.resp(ii).vidSniffs);
        [sniff_frames, sfi,~] = intersect(sniff_frames, fr);
        pfreq = exp.resp(ii).sniffFreq(exp.camTrig(ii).frameInds(sniff_frames));
        traildist = exp.vids(ii).orthogonalDistFromTrail(sniff_frames, 1);
        pos = pos(sfi, :);
        all_pos = cat(1,pos, all_pos);
        all_freq = cat(1, pfreq, all_freq);
        all_dist = cat(1, traildist, all_dist);
        [cm cinds] = getIndexedColors('jet', pfreq);
        for jj=1:length(pfreq)
            line('Parent',gca,'Xdata',pos(jj,1),'Ydata',pos(jj,2),'Marker','.','MarkerSize',8,'Color',cm(cinds(jj),:));
        end
        % to get the scale
        figure; colormap(cm);
        pcolor([pfreq, pfreq]);
        colorbar;
        
        figure; 
        hist(all_freq, 50);
        xlabel('Sniffing Frequency');
    end
end    

