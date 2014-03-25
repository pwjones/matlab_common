function rates = compareOnOffTrailSniffing(exp, trials, threshDist)
% function compareOnOffTrailSniffing(exp, trials)
%
% Compare sniff rates on the trail versus off the trail
% Inputs: experiment structure EXP and the trial numbers to be used
%
% I'm honestly not sure that I have a hypothesis here, just want to see
% if there's a difference between when he's engaged and searching the 
% trail versus when he isn't.

for ii = trials
   fseg = exp.vids(ii).getFollowingSegments([],1,threshDist);
   onframes = [];
   for jj = 1:size(fseg,1)
        onframes = cat(2, onframes, fseg(jj,1):fseg(jj,2));
   end
   offframes = setdiff(1:exp.vids(ii).nFrames, onframes);
   onrate = getSniffRates(exp,ii, onframes, 0);
   offrate = getSniffRates(exp,ii,offframes, 0);
   rates(ii,:) = [onrate, offrate];
end



