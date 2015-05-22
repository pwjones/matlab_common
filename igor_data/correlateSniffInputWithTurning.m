%function correlateSniffInputWithTurning
%
% For each sniff during following, we have the data for the previous 3 sniffs
% in terms of distance from trail, heading, and the sniff to sniff
% differences between those measures.  What we want to do is combine into a
% statistic model, how well can we understand the variation in turning on
% the last intervel from the position of the sniffs before and the change
% in position

dHeading = zeros(size(sniffData.dist_diffs,1), 1);
[paraHeading, absHeading] = convertMouseHeadingAngles(sniffData.mouseHeadingFromOrtho, sniffData.dirToTrail, sniffData.sniffPos);
for ii = 1:1
    dHeading(:,ii) = circ_dist(absHeading(:,ii+3), absHeading(:,ii+1)); 
end
%dHeading = abs(dHeading);
turns = abs(dHeading) > 1;
%turns = turns(:,3);

figure;
subplot(2,3,1);
plot(sniffData.sniffPos(turns,1), dHeading(turns), 'k.');
xlabel('-2 Sniff Position (mm)'); ylabel('Heading Change');
xlim([-15 15]);
subplot(2,3,2);
plot(sniffData.sniffPos(turns,2), dHeading(turns), 'k.');
xlabel('-1 Sniff Position (mm)'); ylabel('Heading Change');
xlim([-15 15]);
subplot(2, 3, 3);
plot(sniffData.sniffPos(turns,3), dHeading(turns), '.k');
xlabel('0 Sniff Position (mm)'); ylabel('Heading Change');
xlim([-15 15]);
subplot(2,3,5);
plot(sniffData.dist_diffs(turns,1), dHeading(turns), 'k.');
xlabel('-2 to -1 \Delta Sniff Position (mm)'); ylabel('Heading Change');
xlim([-10 10]);
subplot(2,3,6);
plot(sniffData.dist_diffs(turns,2), dHeading(turns), 'k.');
xlabel('-1 to 0 \Delta Sniff Position (mm)'); ylabel('Heading Change');
xlim([-10 10]);

posEdges = [-15:-1 1:15];
ddEdges = -10:10;
meanHeading_pos = zeros(length(posEdges),2);
meanHeading_dpos = zeros(length(ddEdges),2);
for jj = 1:2
    [N,bin] = histc(sniffData.sniffPos(turns,jj+1), posEdges);
    for ii = 1:length(posEdges)        
        temp = dHeading(turns, 1);
        temp = temp(bin == ii);
        meanHeading_pos(ii,jj) = nanmean(temp);
    end
    [N,bin] = histc(sniffData.dist_diffs(turns,jj), ddEdges);
    for ii = 1:length(ddEdges)
        temp = dHeading(turns, 1);
        temp = temp(bin == ii);
        meanHeading_dpos(ii,jj) = nanmean(temp);
    end
end

subplot(2,3, 2); ylim([-1 2]);
hold on; plot(posEdges, meanHeading_pos(:,1), 'r', 'LineWidth', 2);
subplot(2, 3, 3); ylim([-1 2]);
hold on; plot(posEdges, meanHeading_pos(:,2), 'r', 'LineWidth', 2);
subplot(2,3,5); ylim([-1 2]);
hold on; plot(ddEdges, meanHeading_dpos(:,1), 'r', 'LineWidth', 2);
subplot(2,3,6); ylim([-1 2]);
hold on; plot(ddEdges, meanHeading_dpos(:,2), 'r', 'LineWidth', 2);