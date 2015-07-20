%% Plotting respiration frequency as a color over the positions
ft = 1:length(exp.resp);
ft = 2;
fr = 600:660;
figure; ah = axes;
for ii = ft
    %fr = 1:exp.vids(ii).nFrames;
    im = ones(exp.vids(ii).height, exp.vids(ii).width);
    im(exp.vids(ii).paths(1).PixelIdxList) = 0;
    im = cat(3, im, ones(exp.vids(ii).height, exp.vids(ii).width), im);
    imshow(im);
    exp.vids(ii).plotPosition(fr, ah, 0, 'k', '.');
    pos = exp.resp(ii).sniffPos;
    sniffFrames = exp.resp(ii).sniffFrames(exp.resp(ii).vidSniffs);
    [sniffFrames, sfi,~] = intersect(sniffFrames, fr);
    sniffFrames = exp.camTrig(ii).frameRange(sniffFrames);
    freq = exp.resp(ii).sniffFreq(exp.camTrig(ii).frameInds(sniffFrames)); 
    pos = pos(sfi, :);
    [cm, cinds] = getIndexedColors('jet', freq, 1);
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