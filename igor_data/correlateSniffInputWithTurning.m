%function correlateSniffInputWithTurning
%
% For each sniff during following, we have the data for the previous 3 sniffs
% in terms of distance from trail, heading, and the sniff to sniff
% differences between those measures.  What we want to do is combine into a
% statistic model, how well can we understand the variation in turning on
% the last intervel from the position of the sniffs before and the change
% in position

resampleSniffData; %run first to randomly resample the behavioral data in
%order to balance mouse positions

posEdges = [-15:-1 1:15];
ddEdges = -9:9;
meanHeading_pos = zeros(length(posEdges),3, size(resamp_idx,3));
meanHeading_dpos = zeros(length(ddEdges),3, size(resamp_idx,3));
meanB = []; meanBint=[];
bturnAll = zeros(size(resamp_idx, 1), size(resamp_idx,3));
for kk=1:size(resamp_idx,3)
    idx = resamp_idx(:,3,kk); %the kk resampling based on most adjacent dd
    
    dHeading = zeros(size(sniffData.dist_diffs(idx,3),1), 1);
    % to avoid discontinuities when animal crosses the trail, convert relative angles
    % to absolute heading angles before taking their distances
    [paraHeading, absHeading] = convertMouseHeadingAngles(sniffData.mouseHeadingFromOrtho(idx,:), sniffData.dirToTrail(idx,:), sniffData.sniffPos(idx,:));
    
    trailDirChange = circ_dist(sniffData.dirToTrail(idx,4), sniffData.dirToTrail(idx,3)); 
    dHeading = circ_dist(absHeading(:,4), absHeading(:,3));
    dHeading = circ_dist(dHeading, trailDirChange); %recorrecting the changes in direction for changes in 
    dHeading = abs(dHeading);
    bturn = logical(sniffData.bturn(idx));
    
    for jj = 1:3
        [N,bin] = histc(sniffData.sniffPos(idx,jj+1), posEdges);
        for ii = 1:length(posEdges)
            temp = dHeading(:, 1);
            temp = temp(bin == ii);
            meanHeading_pos(ii,jj, kk) = nanmean(temp);
        end
        [N,bin] = histc(sniffData.dist_diffs(idx,jj), ddEdges);
        for ii = 1:length(ddEdges)
            temp = dHeading(:, 1);
            temp = temp(bin == ii);
            meanHeading_dpos(ii,jj, kk) = nanmean(temp);
        end
    end
    
    X = [sniffData.dist_diffs(idx,3), sniffData.dist_diffs(idx,2), sniffData.sniffPos(idx, 4), sniffData.sniffPos(idx, 3) ones(length(idx),1)];
    [b,bint,r,rint,stats] = regress(dHeading,X);
    meanB = cat(2, meanB, b);
    meanBint = cat(3, meanBint, bint);
end
meanHeading_pos = nanmean(meanHeading_pos, 3);
meanHeading_dpos = nanmean(meanHeading_dpos, 3);

figure;
ah = [];
subplot(2,3,1);
ah(1) =plot(sniffData.sniffPos(idx,2), dHeading(:), '.k'); hold on;
plot(sniffData.sniffPos(idx(bturn),2), dHeading(bturn), '.b');
xlabel('-2 Sniff Position (mm)'); ylabel('Heading Change');
xlim([-15 15]);
subplot(2,3,2);
ah(2) = plot(sniffData.sniffPos(idx,3), dHeading(:), '.k'); hold on;
plot(sniffData.sniffPos(idx(bturn),3), dHeading(bturn), '.b');
xlabel('-1 Sniff Position (mm)'); ylabel('Heading Change');
xlim([-15 15]);
subplot(2, 3, 3);
ah(3) = plot(sniffData.sniffPos(idx,4), dHeading(:), '.k'); hold on;
xlabel('0 Sniff Position (mm)'); ylabel('Heading Change');
xlim([-15 15]);
subplot(2,3,5);
ah(5) = plot(sniffData.dist_diffs(idx,2), dHeading(:), 'k.');
xlabel('-2 to -1 \Delta Sniff Position (mm)'); ylabel('Heading Change');
xlim([-10 10]);
subplot(2,3,6);
ah(6) = plot(sniffData.dist_diffs(idx,3), dHeading(:), 'k.'); hold on;
plot(sniffData.dist_diffs(idx(bturn),3), dHeading(bturn), '.b');
xlabel('-1 to 0 \Delta Sniff Position (mm)'); ylabel('Heading Change');
xlim([-10 10]);

subplot(2,3, 2); ylim([0 2]);
hold on; plot(posEdges, meanHeading_pos(:,1), 'r', 'LineWidth', 2);
subplot(2, 3, 3); ylim([0 2]);
hold on; plot(posEdges, meanHeading_pos(:,2), 'r', 'LineWidth', 2);
subplot(2,3,5); ylim([0 2]);
hold on; plot(ddEdges, meanHeading_dpos(:,1), 'r', 'LineWidth', 2);
subplot(2,3,6); ylim([0 2]);
hold on; plot(ddEdges, meanHeading_dpos(:,2), 'r', 'LineWidth', 2);