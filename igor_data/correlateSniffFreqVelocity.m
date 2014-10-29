function [SF, NV, NV_filt, NA] = correlateSniffFreqVelocity(exp_struct)
%% Look at the correlation between velocity (body and nose) with sniff frequency

% analysis parameters
mm_conv = 1.16; %mm/px linear
followDistThresh = 20; %px

all_freq = []; all_noseVel_filt = []; all_noseVel = [];
mov_thresh = 50;
for ii = 1:length(exp_struct.resp)
    % select the frames to analyze - when the animal is following the trail
    sniffFrames = exp_struct.resp(ii).sniffFrames(exp_struct.resp(ii).vidSniffs);
    sniffDists = exp_struct.vids(ii).orthogonalDistFromTrail(sniffFrames, 1);
    followi = sniffDists <= followDistThresh & sniffDists >= -followDistThresh;
    sniffFrames = sniffFrames(followi); 
    noseVel = exp_struct.vids(ii).noseVel; %get the nose/body velocities
    bodyVel = exp_struct.vids(ii).bodyVel(:,1);
    noseVel_filt = gaussianFilter(noseVel, 1, 'conv'); %smoother versions - vels tend to look messy
    bodyVel_filt = gaussianFilter(bodyVel, 1, 'conv'); 
    noseVel = noseVel(sniffFrames); noseVel_filt = noseVel_filt(sniffFrames);%select the frames
    bodyVel = bodyVel(sniffFrames); bodyVel_filt = bodyVel_filt(sniffFrames);
    sniffFrames = exp_struct.camTrig(ii).frameRange(sniffFrames); %index frames collected rather than those analysed
    freq = exp_struct.resp(ii).sniffFreq(exp_struct.camTrig(ii).frameInds(sniffFrames));
    norm_freq = freq ./ nanmax(freq);
    
    % unit conversions - from px/frame to mm/sec
    fps = exp_struct.vids(ii).frameRate;
    noseVel = mm_conv * fps * noseVel; noseVel_filt = mm_conv * fps * noseVel_filt; 
    bodyVel = mm_conv * fps * bodyVel; bodyVel_filt = mm_conv * fps * bodyVel_filt; 
    
    % Let's try to see what things look like excluding the really low vel points
    moving = noseVel_filt >= mov_thresh;
    freq = freq(moving); noseVel = noseVel(moving); bodyVel = bodyVel(moving);
    noseVel_filt = noseVel_filt(moving);
    
    % save them
    all_freq = cat(1, all_freq, freq(:));
    all_noseVel_filt = cat(1,all_noseVel_filt, noseVel_filt(:));
    all_noseVel = cat(1,all_noseVel, noseVel(:));    
end

SF = all_freq;
NV = all_noseVel;
NV_filt = all_noseVel_filt;
NA = [0; diff(NV_filt)];
