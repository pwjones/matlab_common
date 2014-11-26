%% Script to look at the distribution of sniff positions relative to the trail.

mm_conv = .862; %mm/px linear
f1= figure;
all_sniffDists = []; all_followDists = [];
thresh_dist = 20;
for ii = 1:length(exp.resp)
    exp.vids(ii).makePathsSkel();
    % Velocities
    noseVel = exp.vids(ii).noseVel * mm_conv * exp.vids(ii).frameRate; %get the nose/body velocities
    nosePos = exp.vids(ii).nosePos ;
    bodyVel = exp.vids(ii).bodyVel(:,1) * mm_conv * exp.vids(ii).frameRate;
    noseVel_filt = gaussianFilter(noseVel, 3, 'conv'); %smoother versions - vels tend to look messy
    bodyVel_filt = gaussianFilter(bodyVel, 3, 'conv'); 
    % select the frames to analyze - when the animal is following the trail
    allDists = exp.vids(ii).orthogonalDistFromTrail(1:exp.vids(ii).nFrames, 1);
    followAll = allDists <= 20 & allDists > -20;
    moving = noseVel_filt >= 50;
    followDists = allDists(followAll & moving);
    all_followDists = cat(1, all_followDists, followDists(:));
    % when the animal is following and sniffing
    sniffFrames = exp.resp(ii).sniffFrames(exp.resp(ii).vidSniffs);
    sniffDists = exp.vids(ii).orthogonalDistFromTrail(sniffFrames, 1);
    followi = sniffDists <= 20 & sniffDists >= -20;
    moving = noseVel_filt(sniffFrames) >= 50;
    sniffFrames = sniffFrames(followi & moving); 
    
    noseVel = noseVel(sniffFrames); noseVel_filt = noseVel_filt(sniffFrames);%select the frames
    bodyVel = bodyVel(sniffFrames); bodyVel_filt = bodyVel_filt(sniffFrames);
    following_sniffDists = sniffDists(followi);
    all_sniffDists = cat(1, all_sniffDists, following_sniffDists(:));
    sniffFrames = exp.camTrig(ii).frameRange(sniffFrames); %index frames collected rather than those analysed
    freq = exp.resp(ii).sniffFreq(exp.camTrig(ii).frameInds(sniffFrames));
end

xbins = linspace(-thresh_dist, thresh_dist, 40);
sniff_hist = histc(all_sniffDists, xbins);
all_hist = histc(all_followDists, xbins);
figure(f1);
plot(xbins, 100*all_hist./sum(all_hist), 'k', 'LineWidth',2); hold on;
plot(xbins, 100*sniff_hist./sum(sniff_hist), 'g', 'LineWidth', 2);
xlabel('Distance from trail (px)');
ylabel('Percent of following time');
legendstr = {'All Frames', 'Sniff Frames'};
legend(legendstr, 'Location', 'NorthEast');
title(extractMouseNameFromFN(exp.vids(1).videoFN));
