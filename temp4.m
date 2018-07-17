trailDirChange = circ_dist(sniffData.dirToTrail(:,4), sniffData.dirToTrail(:,3)); 
turningTrailRelative = circ_dist(paraHeading, trailDirChange);
figure;
for ii=1:200
    %polarplot(sniffData.mouseHeadingFromOrtho(ii,:), sniffData.sniffPos(ii,:));
    %polarplot(paraHeading(ii,:), abs(sniffData.sniffPos(ii,:)));
    ind = (ii-1)*4 + (1:4);
    polarplot(turningTrailRelative(ii,:)', ind);
    hold on;
end

v=perMouseData(1).nose_vel{178};
nv=[];
for jj=1:length(v)
    nv = cat(1, nv, v{jj});
end