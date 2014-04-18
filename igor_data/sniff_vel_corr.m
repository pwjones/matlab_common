%% Look at the correlation between velocity (body and nose) with sniff frequency
mm_conv = 1.16; %mm/px linear
f1= figure;
f2 = figure;
f3 = figure;
all_freq = []; all_noseVel_filt = []; all_noseVel = [];
mov_thresh = 50;
for ii = 1:length(exp.resp)
    % select the frames to analyze - when the animal is following the trail
    sniffFrames = exp.resp(ii).sniffFrames(exp.resp(ii).vidSniffs);
    sniffDists = exp.vids(ii).orthogonalDistFromTrail(sniffFrames, 1);
    followi = sniffDists <= 20 & sniffDists >= -20;
    sniffFrames = sniffFrames(followi); 
    noseVel = exp.vids(ii).noseVel; %get the nose/body velocities
    bodyVel = exp.vids(ii).bodyVel(:,1);
    noseVel_filt = gaussianFilter(noseVel, 3, 'conv'); %smoother versions - vels tend to look messy
    bodyVel_filt = gaussianFilter(bodyVel, 3, 'conv'); 
    noseVel = noseVel(sniffFrames); noseVel_filt = noseVel_filt(sniffFrames);%select the frames
    bodyVel = bodyVel(sniffFrames); bodyVel_filt = bodyVel_filt(sniffFrames);
    sniffFrames = exp.camTrig(ii).frameRange(sniffFrames); %index frames collected rather than those analysed
    freq = exp.resp(ii).sniffFreq(exp.camTrig(ii).frameInds(sniffFrames));
    %noseVel_filt = noseVel_filt ./ nanmax(noseVel_filt); %normalize
    %bodyVel_filt = bodyVel_filt ./ nanmax(bodyVel_filt); 
    norm_freq = freq ./ nanmax(freq);
    % unit conversions - from px/frame to mm/sec
    fps = exp.vids(ii).frameRate;
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
    %correlation/fit
    %[p, S] = polyfit(noseVel, freq, 1);
    %fitx = [0 max(noseVel)];
    %fity = polyval(p, fitx);
    %r = corr(noseVel, freq);
    figure(f1);
    subplot (3, ceil(length(exp.vids)/3), ii);
    plot(noseVel_filt, freq, 'b.', 'MarkerSize', 12);
    %plot(noseVel, freq, 'g.', 'MarkerSize', 12);
    %hold on; plot(bodyVel, freq,'mx');
    %text(50, 8, ['y = ' num2str(p(1)) 'x + ' num2str(p(2)) ', r = ' num2str(r)]); 
    %plot(fitx, fity, 'b', 'LineWidth', 1);
    xlim([50 400]); ylim([5 18]);
    xlabel('Vel (mm/sec)'); ylabel('Sniff freq (Hz)');
    lg_str = {'Nose Vel', 'Body Vel'};
    %legend(lg_str, 'Location', 'SouthEast');
    
    %figure(f2);
    %subplot(3, ceil(length(exp.vids)/3), ii);
    %plot(1:15, 1:15,'k--'); hold on;
    %plot(bodyVel, noseVel, 'r.'); xlabel('Body Vel'); ylabel('Nose Vel');
    %xlim([0 400]); ylim([0 400]);
     
    figure(f3);
    subplot(3, ceil(length(exp.vids)/3), ii);
    x = 1:length(freq);
    %plot(x, noseVel_filt, 'b-', x, bodyVel_filt, 'm-', x, 50*norm_freq, 'r');
    lh = plot(x, noseVel_filt, 'b-', x, freq*3, 'r');
    set(lh, 'LineWidth', 1);
    xlabel('Sniff Number during following'); ylabel('Velocity(mm/s) and Freq (Hz * 3)');
end

f4 = figure;
plot(all_noseVel_filt, all_freq, 'b.', 'MarkerSize', 12);
xlim([50 400]); ylim([5 18]);
xlabel('Vel (mm/sec)'); ylabel('Sniff freq (Hz)');
%correlation/fit
[p, S] = polyfit(all_noseVel_filt, all_freq, 1)
fitx = [mov_thresh max(all_noseVel_filt)];
fity = polyval(p, fitx);
r = corr(all_noseVel, all_freq);
hold on; plot(fitx, fity, 'r');
% Build a binned curve
bs = 25;
vel_bin = mov_thresh:bs:400;
binned_freq = zeros(length(vel_bin)-1,1);
for jj = 1:(length(vel_bin)-1)
    bi = all_noseVel_filt >= vel_bin(jj) & all_noseVel_filt < vel_bin(jj+1);
    binned_freq(jj) = nanmean(all_freq(bi));
end
vel_x = vel_bin(1:end-1) + bs/2;
hold on; plot(vel_x, binned_freq, 'r', 'LineWidth',1);