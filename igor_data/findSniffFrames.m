function exp = findSniffFrames(exp)
% function exp = findSniffFrames(exp)
%
% processing the integrated sniffing and video
% 


% Let's get the frames when the animal sniffed
for ii=1:length(exp.resp)
    % These are the frames that have been collected - not necessarily
    % analyzed in videos
    frame_t = exp.camTrig(ii).time(exp.camTrig(ii).frameInds);
    sniff_times = exp.resp(ii).time(exp.resp(ii).sniffVect);
    sniff_frames = zeros(length(sniff_times),1);
    for jj=1:length(sniff_times)
        sf_temp = find(sniff_times(jj) >= frame_t, 1, 'last');
        if isempty(sf_temp)
            sniff_frames(jj) = 0;
        else
            sniff_frames(jj) = sf_temp;
        end
    end
    exp.resp(ii).sniffTimes = sniff_times;
    %adjust the frame for when the analysis starts
    exp.resp(ii).sniffFrames = sniff_frames - double(exp.camTrig(ii).vidStartFrame) + 1; 
    exp.resp(ii).vidSniffs = find(exp.resp(ii).sniffFrames >= 1 & exp.resp(ii).sniffFrames <= exp.vids(ii).nFrames);
    vsi = exp.resp(ii).vidSniffs; %the indices within the total sniffs, of the ones that we have analyzed video for
    exp.resp(ii).nSniffDuringVid = length(exp.resp(ii).vidSniffs); 
    exp.resp(ii).sniffPos = NaN * zeros(exp.resp(ii).nSniffDuringVid,2); 
    exp.resp(ii).sniffPos = exp.vids(ii).nosePos(exp.resp(ii).sniffFrames(vsi),:);
    exp.resp(ii).vidTime = exp.resp(ii).time - exp.resp(ii).time(exp.camTrig(ii).analyzedFrameInds(1));
    exp.resp(ii).vidSniffTimes = exp.resp(ii).sniffTimes - exp.resp(ii).time(exp.camTrig(ii).analyzedFrameInds(1));
end
