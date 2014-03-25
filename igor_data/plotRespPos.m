%% Plotting respiration frequency as a color over the positions
ft = 1:length(exp.vids);
ft = 1;
fr = [];
fr = 1100:1580;
for ii = ft
    %fr = 1:exp.vids(ii).nFrames;
    exp.vids(ii).plotPosition(fr, [], 0, 'k', '.');
    pos = exp.resp(ii).sniffPos;
    sniff_frames = exp.resp(ii).sniffFrames(exp.resp(ii).vidSniffs);
    [sniff_frames, sfi,~] = intersect(sniff_frames, fr);
    freq = exp.resp(ii).sniffFreq(exp.camTrig(ii).frameInds(sniff_frames)); 
    pos = pos(sfi, :);
    [cm cinds] = getIndexedColors('jet', freq);
%     minfreq = min(freq)
%     maxfreq = max(freq)
    for jj=1:length(freq)
        line('Parent',gca,'Xdata',pos(jj,1),'Ydata',pos(jj,2),'Marker','.','MarkerSize',8,'Color',cm(cinds(jj),:));
    end
    % to get the scale
    figure; colormap(cm);
    pcolor([freq, freq]);
    colorbar;
end