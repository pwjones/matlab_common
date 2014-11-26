% Let's get the frames when the animal sniffed
for ii=1:length(exp.resp)
    frame_t = exp.camTrig(ii).time(exp.camTrig(ii).frameInds);
    sniff_times = exp.resp(ii).time(exp.resp(ii).sniffVect);
    sniffFrames = zeros(length(sniff_times),1)*NaN;
    for jj=1:length(sniff_times)
        tempf =  find(frame_t >= sniff_times(jj), 1, 'first');
        if ~isempty(tempf)
            sniffFrames(jj) = tempf;
        else 
            sniffFrames(jj) = length(frame_t);
        end
    end
    exp.resp(ii).sniff_times = sniff_times;
    exp.resp(ii).sniff_frames = sniffFrames;
    exp.resp(ii).sniff_pos = exp.vids(ii).nosePos(sniffFrames,:);
end

%% Plot the mouse positions
vidi = 1:length(vids);
vidi = 1:25;
%vidi = 24:28;
for ii = 1:length(vidi)
    %exp.vids(ii).plotNosePosition([]);
    vids(ii).plotNosePosition([]);
    %vids(ii).plotFollowing([],15,'');
    %exp.vids(ii).plotPosition([exp.resp(ii).sniffFrames(exp.resp(ii).vidSniffs)], ah, 0, 'b', '.');
end
%% Just print the video file names
for ii = 1:length(vids)
    disp(vids(ii).videoFN);
    %exp.vids(ii).blobID(100,:)
    %exp.vids(ii).fcPeriod = 60;
    vids(ii).save;
    
end
%% Clear the tracked points that are actually LED
for ii = 1:length(vids)
    np = vids(ii).nosePos;
    corner = np(:,1) >= 1206 & np(:, 2) >= 954;
    vids(ii).nosePos(np, :) = [NaN NaN];
    vids(ii).computeVelocity([]);
    vids(ii).save;
end

%% Look at the correlation between velocity (body and nose) with sniff frequency
f1= figure;
f2 = figure;
for ii = 1:length(exp.resp)
    % select the frames to analyze
    sniffFrames = exp.resp(ii).sniffFrames(exp.resp(ii).vidSniffs);
    sniffDists = exp.vids(ii).orthogonalDistFromTrail(sniffFrames, 1);
    followi = sniffDists <= 20;
    sniffFrames = sniffFrames(followi); 
    %select values
    noseVel = exp.vids(ii).noseVel(sniffFrames);
    bodyVel = exp.vids(ii).bodyVel(sniffFrames);
    bodyVel = sqrt(sum(bodyVel.^2,2)); %make it just a magnitude
    sniffFrames = exp.camTrig(ii).frameRange(sniffFrames);
    freq = exp.resp(ii).sniffFreq(exp.camTrig(ii).frameInds(sniffFrames));
    size(noseVel)
    size(freq)
    figure(f1);
    subplot (3, ceil(length(exp.vids)/3), ii);
    plot(noseVel, freq, 'bo');
    hold on; plot(bodyVel, freq,'gx');
    xlim([0 7]);
    figure(f2);
    subplot(3, ceil(length(exp.vids)/3), ii);
    plot(1:15, 1:15,'k--'); hold on;
    plot(bodyVel, noseVel, 'r.'); xlabel('Body Vel'); ylabel('Nose Vel');
    xlim([0 7]); ylim([0 7])
end

%% Plotting respiration frequency as a color over the positions
ft = 1:length(exp.resp);
ft = 3;
fr = 1:1500;
for ii = ft
    %fr = 1:exp.vids(ii).nFrames;
    exp.vids(ii).plotPosition(fr, [], 0, 'k', '.');
    pos = exp.resp(ii).sniffPos;
    sniffFrames = exp.resp(ii).sniffFrames(exp.resp(ii).vidSniffs);
    [sniffFrames, sfi,~] = intersect(sniffFrames, fr);
    sniffFrames = exp.camTrig(ii).frameRange(sniffFrames);
    freq = exp.resp(ii).sniffFreq(exp.camTrig(ii).frameInds(sniffFrames)); 
    pos = pos(sfi, :);
    [cm cinds] = getIndexedColors('jet', freq, 1);
    %minfreq = min(freq)
    %maxfreq = max(freq)
    for jj=1:length(freq)
        line('Parent',gca,'Xdata',pos(jj,1),'Ydata',pos(jj,2),'Marker','.','MarkerSize',8,'Color',cm(cinds(jj),:));
    end
    % to get the scale
    figure; colormap(cm);
    pcolor([freq, freq]);
    colorbar;
end


%%
for ii = 1:length(exp.vids)
    %exp.vids(ii).detectPaths(1,1,1);
    %exp.vids(ii).refinePaths(1);
    %exp.vids(ii).refinePaths(2);
    exp.vids(ii).save;
end

%%
disp('Checking for missing frames in videos');
for ii = 1:length(exp.vids)
    missing = exp.vids(ii).isMissingFrame();
    if missing
        disp(sprintf('%s is missing a frame.', exp.vids(ii).videoFN));
    end
end
%% Let's try to look at the following and sniffing together.

  ft = [3 4 6];
  %ft = 7;
  %fr = 330:900;
  for ii = ft
     fr = 1:exp.vids(ii).nFrames;
     exp.vids(ii).plotFollowing(fr, 15,'');
     fr = intersect(fr, exp.resp(ii).sniffFrames(exp.resp(ii).vidSniffs));
     exp.vids(ii).plotNosePosition(fr, gca, 0, 'y', '.');
  end
  
 %% just get some paths

 vidi = 1:length(exp.vids);
 for ii = vidi
        exp.vids(ii).detectRefinePaths(1,1,1);
        exp.vids(ii).save();
 end
 
 %% just get some paths
 vidi = length(vids):-1:1;
 vidi = 1:length(vids);
 %vidi = 22:-1:1;
 %vidi = 1:12;
 %vidi = 9;
 for ii = vidi
        vids(ii).detectRefinePaths(1,1,1);
        vids(ii).save();
 end
 
 %%
 nose_mean = NaN*ones(size(nose_px,1),1);
 for ii = 1:size(nose_px,1)
    nose_mean(ii) = nanmean(nose_px(ii,:));
 end
 
 %% Plot jposition/velocity for videos
 for ii =1:length(vids)
     vids(ii).plotVelocity([]);
 end
 
 %% Plot the nose velocity as traces
 vidi = 1:length(vids);
 figure; hold on;
 for ii = vidi
     vel = sqrt(sum(vids(ii).noseVel.^2,2));
     t = vids(ii).times(:);
     plot(t,vel, 'b');
 end
 
 %% Compute velocity for a set of videos
 vidi = 1:length(vids);
 for ii=vidi
    vids(ii).computeVelocity(1:vids(ii).nFrames); 
 end
 
 %% Looking at the position dependence of the nose brightness
 irange = 1:mt.nFrames;
 ft = 3;
 figure; imshow(frame); colormap('gray');
 pos = mt.nosePos(irange,:);
 plot_fr = 1:2600;
 [cm cinds] = getIndexedColors('jet', meanlum(plot_fr));
 for ii=1:length(plot_fr)
     jj = plot_fr(ii);
     if ~isnan(meanlum(jj))
        line('Parent',gca,'Xdata',pos(jj,1),'Ydata',pos(jj,2),'Marker','.','MarkerSize',8,'Color',cm(cinds(ii),:));
     end
 end
 % to get the scale
 figure; colormap(cm);
 pcolor([meanlum, meanlum]);
 colorbar;
 
 
 %% Plot the velocities of the body and nose                                                                                                                                                                                                                                                                                                                                                                                                                          
figure; 
for ii = 1:length(vids)
    subplot (3, ceil(length(vids)/3), ii);
    bc = vids(ii).bodyCOM;
    np = vids(ii).nosePos;
    lh = plot(bc(:,1), bc(:,2)); 
    hold on; lh2 = plot(np(:,1), np(:,2), 'Color', [0 .7 0]); 
end

figure; 
for ii = 1:length(vids)
    subplot (3, ceil(length(vids)/3), ii);
    %np = vids(ii).findNose(1:vids(ii).nFrames);
    np = vids(ii).nosePos;
    nf = vids(ii).nFrames;
    disp('\n');
    disp(['Number of frames without nose: ' num2str(sum(isnan(np(:,1))))]);
    vids(ii).computeVelocity(1:vids(ii).nFrames);
    bv = sqrt(sum(vids(ii).bodyVel .^2, 2));
    nv = sqrt(sum(vids(ii).noseVel.^2, 2));
    disp(['Number of frames without nose velocity: ' num2str(sum(isnan(nv)))]);
    lh = plot(vids(ii).times(:), bv, vids(ii).times(:), nv); 
    set(lh(2), 'LineWidth', 1);%plot(vids(ii).times(:), bv);
    disp(['Total frames: ' num2str(nf)]);
    disp(['Percent without nose: ' num2str(sum(isnan(np(:,1)))./nf*100)]);
    disp(['Percent without nose velocity: ' num2str(sum(isnan(nv)./nf*100))]);
    %vids(ii).fcPeriod = 60;
    %vids(ii).isMissingFrame();
    %vids(ii).plotNosePosition([]);
end
 