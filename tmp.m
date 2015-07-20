%% Plotting respiration frequency as a color over the positions
ft = 1:length(exp.resp);
ft = 2;
fr = [];
plotblue = [.3, .6, 1];
for ii = ft
    %fr = 1:exp.vids(ii).nFrames;
    fr= 500:1000;
    bgIm = exp.vids(ii).plotPaths(); 
    figure; imshow(255*bgIm); hold on;
    np = exp.vids(ii).nosePos;
    hold on; plot(np(:,1), np(:,2), '-', 'Marker', '.','MarkerSize', 24, 'Color', plotblue, 'LineWidth', 3);
    
    pos = exp.resp(ii).sniffPos;
    sniffFrames = exp.resp(ii).sniffFrames(exp.resp(ii).vidSniffs);
    [sniffFrames, sfi,~] = intersect(sniffFrames, fr);
    sniffFrames = exp.camTrig(ii).frameRange(sniffFrames);
    freq = exp.resp(ii).sniffFreq(exp.camTrig(ii).frameInds(sniffFrames)); 
    pos = pos(sfi, :);
    %[cm, cinds] = getIndexedColors('jet', freq, 1);
    %minfreq = min(freq)
    %maxfreq = max(freq)
    for jj=1:length(freq)
        line('Parent',gca,'Xdata',pos(jj,1),'Ydata',pos(jj,2),'Marker','.','MarkerSize',26,'Color','r');
    end
    turni = exp.vids(ii).findFollowingTurns([], 1, 20, -3:10);
    tp = np(turni,:);
    plot(tp(:,1), tp(:,2), 'y.', 'MarkerSize', 20);
    % to get the scale
    %figure; colormap(cm);
    %pcolor([freq, freq]);
    %colorbar;
end