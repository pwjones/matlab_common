%function correlateSniffInputWithTurning2
%
% This version only looks at actual detected turns!
%
% For each sniff during following, we have the data for the previous 3 sniffs
% in terms of distance from trail, heading, and the sniff to sniff
% differences between those measures.  What we want to do is combine into a
% statistic model, how well can we understand the variation in turning on
% the last intervel from the position of the sniffs before and the change
% in position

[~,pre] = convertMouseHeadingAngles(sniffData.turnTrig_preTurnHeadings, sniffData.turnTrig_preTurnDirToTrail, sniffData.turnTrig_preTurnPos); 
[~,post] = convertMouseHeadingAngles(sniffData.turnTrig_postTurnHeadings, sniffData.turnTrig_postTurnDirToTrail, sniffData.turnTrig_postTurnPos); 
dHeading = circ_dist(post, pre); 
%dHeading = dHeading;

figure;
subplot(2,2, 1);
plot(sniffData.turnTrig_sniffPos(:,end-1), dHeading, 'k.');
xlabel('-1 Sniff Position (mm)'); ylabel('Heading Change');
xlim([-15 15]); %ylim([1 2]);
subplot(2, 2, 2);
plot(sniffData.turnTrig_sniffPos(:,end), dHeading, 'k.');
xlabel('0 Sniff Position (mm)'); ylabel('Heading Change');
xlim([-15 15]); %ylim([1 2]);
subplot(2,2,3);
plot(sniffData.turnTrig_preTurnDistDiff(:,end-1), dHeading, 'k.');
xlabel('-2 to -1 \Delta Sniff Position (mm)'); ylabel('Heading Change');
xlim([-10 10]); %ylim([1 2]);
subplot(2,2,4);
plot(sniffData.turnTrig_preTurnDistDiff(:,end), dHeading, 'k.');
xlabel('-1 to 0 \Delta Sniff Position (mm)'); ylabel('Heading Change');
xlim([-10 10]); %ylim([1 2]);

posEdges = -20:20;
ddEdges = -15:15;
for jj = 1:2
    for ii = 1:length(posEdges)        
        [N,bin] = histc(sniffData.turnTrig_sniffPos(:,jj+2), posEdges);
        temp = dHeading(bin == ii);
        meanHeading_pos(ii,jj) = nanmean(temp);
    end
    for ii = 1:length(ddEdges)
        [N,bin] = histc(sniffData.turnTrig_preTurnDistDiff(:,jj+1), ddEdges);
        temp = dHeading(bin == ii);
        meanHeading_dpos(ii,jj) = nanmean(temp);
    end
end

subplot(2,2, 1);
hold on; plot(posEdges, meanHeading_pos(:,1), 'r', 'LineWidth', 2);
xlim([-15 15]); 
subplot(2, 2, 2);
hold on; plot(posEdges, meanHeading_pos(:,2), 'r', 'LineWidth', 2);
xlim([-15 15]);
subplot(2,2,3);
hold on; plot(ddEdges, meanHeading_dpos(:,1), 'r', 'LineWidth', 2);
xlim([-10 10]);
subplot(2,2,4);
hold on; plot(ddEdges, meanHeading_dpos(:,2), 'r', 'LineWidth', 2);
xlim([-10 10]);