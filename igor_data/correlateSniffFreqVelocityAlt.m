function [SF, NV, NV_filt, NA, SF_ISI, cc] = correlateSniffFreqVelocityAlt(exp_struct, followDistThresh)

% analysis parameters
%mm_conv = 1.16; %mm/px linear
mm_conv = .862; %mm/px linear
%followDistThresh = 20; %px
pb = 0; 

all_freq = []; all_noseVel_filt = []; all_noseVel = []; all_freq_isi = [];
mov_thresh = 25;
for ii = 1:length(exp_struct.resp)
    % In the other way of analyzing this data, I find every frame where there is a sniff
    % and then look at those frames specifically.  I'm not convinced that I should only 
    % look at the frames where there is a sniff, so I'm relaxing that and using every frame.
    % select the frames to analyze - when the animal is following the trail
    
    %sniffFrames = exp_struct.resp(ii).sniffFrames(exp_struct.resp(ii).vidSniffs);
    sniffFrames = 1:min(exp_struct.vids(ii).nFrames, exp_struct.camTrig(ii).nFrames);
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
    ISI_offset = int32((exp_struct.resp(ii).ISIrate_t(1) - exp_struct.resp(ii).ts) /exp_struct.resp(ii).dt);
    freq_isi = [0; exp_struct.resp(ii).ISIrate(exp_struct.camTrig(ii).frameInds(sniffFrames(1:end-1))-double(ISI_offset))];
    %norm_freq = freq ./ nanmax(freq);
    
    % unit conversions - from px/frame to mm/sec
    fps = exp_struct.vids(ii).frameRate;
    noseVel = mm_conv * fps * noseVel; noseVel_filt = mm_conv * fps * noseVel_filt; 
    bodyVel = mm_conv * fps * bodyVel; bodyVel_filt = mm_conv * fps * bodyVel_filt; 
    
    % Let's try to see what things look like excluding the really low vel points
    moving = noseVel_filt >= mov_thresh;
    freq = freq(moving); noseVel = noseVel(moving); bodyVel = bodyVel(moving); freq_isi = freq_isi(moving);
    noseVel_filt = noseVel_filt(moving);
    
    % save them
    all_freq = cat(1, all_freq, freq(:));
    all_noseVel_filt = cat(1,all_noseVel_filt, noseVel_filt(:));
    all_noseVel = cat(1,all_noseVel, noseVel(:));    
    all_freq_isi = cat(1, all_freq_isi, freq_isi(:));
end

SF = all_freq;
NV = all_noseVel;
NV_filt = all_noseVel;
NA = [0; diff(NV_filt)];
SF_ISI = all_freq_isi;
cc = [];

%Compute the cross correlation of the time series as in the wikipedia article on cross correlation 
xwind = 50;
[c, lags] = xcorr(NV_filt-nanmean(NV_filt), SF-nanmean(SF), xwind, 'biased');
c = c./(nanstd(NV_filt)*nanstd(SF));
cc.corr = c;
cc.lags = lags;

if pb
    figure;
    stem(lags, c);
end

% Shuffle correction for cross correlation 
%     n_shuffle = 200;
%     c_shuffle = zeros(length(lags), n_shuffle);
%     for jj = 1:n_shuffle
%         sel = randsample(length(SF), length(SF));
%         c_shuffle(:,jj) = xcorr(NV_filt, SF(sel), xwind, 'coeff');
%     end
%     c_shuffle = mean(c_shuffle,2);
    %stem(lags, c-c_shuffle);
    

    