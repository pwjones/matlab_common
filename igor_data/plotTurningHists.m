function plotTurningHists(sniffData)
%function plotTurningRoseHists(sniffData)
%
% Makes a few plots looking at the pre and 
% post turn angles and turn magnitudes relative to 
% the trail and the magnitude of the distance change 
% pre-turn.

%                   dist_diffs: [36547x1 double]
%                   bturn:
%                turnDirections: [36547x1 double]
%                  mouseHeading: [36547x1 double]
%      turnTrig_preTurnHeadings: [9204x1 double]
%     turnTrig_turnTrajectories: [21x9204 double]
%             turnTrig_turnDirs: [9204x1 double]
%     turnTrig_postTurnHeadings: [9204x1 double]
%      turnTrig_preTurnDistDiff: [9204x1 double]
%              turnTrig_turnLag: []
%        turnTrig_sniffHeadings: [9204x1 double]
       
       
% Figure 1 - Distance changes, turning probabilities
% ---------------------------------------------------
turnDirections = sniffData.turnDirections./abs(sniffData.turnDirections); %convert to signed 1's.
edges = -9.5:9.5;
for jj = 1:size(sniffData.dist_diffs,2)
    [N,bin] = histc(sniffData.dist_diffs(:,jj), edges);
    nbins = length(N);
    turn_counts = zeros(1,nbins);
    turnTowards= zeros(1,nbins);
    turnTowardsProp = zeros(1,nbins);
    for ii = 1:nbins
        turn_counts(ii) = sum(sniffData.bturn(bin == ii));
        temp = turnDirections(bin == ii); % turns after each concentration change
        turnTowards(ii) = sum(temp(temp == 1)); %The number of turns towards the trail
        turnTowardsProp(ii) = turnTowards(ii)/turn_counts(ii);
        %headings = mouseHeading(bin == ii);
        %meanHeadings(ii) = mean(headings);
        %stdHeadings(ii) = std(headings);
    end
    turnp(:,jj) = turn_counts(:)./N(:);
end

turnp(isnan(turnp)) = 0; %replace possible NaNs

% similar, but doing it using the turning triggered data.
for ii = 1:nbins
    for jj = 1:3
        [N_s, bin_s] = histc(sniffData.turnTrig_preTurnDistDiff(:,jj), edges);
        turns = sniffData.turnTrig_turnDirs(bin_s==ii);
        turnToP(ii,jj) = sum(turns>0)./numel(turns);
    end
end
 
figure;
subplot(1,3,1);
%bar(edges(1:end-1), N, 'histc');
stairs(edges, N);
xlabel('\Delta Distance from Trail (mm)');
ylabel('Counts');
xlim([-10 10]);

subplot(1,3,2);
%bar(edges(1:end-1),turnp, 'histc');
stairs(edges',turnp);
xlabel('\Delta Distance from Trail (mm)');
ylabel('Turning Probability');
xlim([-10 10]);
ylim([0 .5]);

subplot(1,3,3);
%bar(edges(1:end-1),turnTowardsProp, 'histc');
%bar(edges, turnTowardsProp, 'histc');
stairs(edges', turnToP);
xlabel('\Delta Distance from Trail (mm)');
ylabel('Prop Turns Towards Trail');
legend({'-3to-4', '-3to-2','-2to-1'});
xlim([-10 10]);

% Figure 2 - Headings Pre and Post Turn
% -------------------------------------
figure;
subplot(2,3,1);
cats = [-10 -2.5 2.5 10];
[Ncat, catbin] = histc(sniffData.turnTrig_preTurnDistDiff(:,3), cats);
pre = sniffData.turnTrig_preTurnHeadings(:,1);
post = sniffData.turnTrig_postTurnHeadings(:,1);
rh = rose2(mod(pre(catbin == 1)+2*pi, 2*pi));
title('Approach - Pre');
%---------------
subplot(2,3,2);
rh = rose2(mod(pre(catbin == 2)+2*pi, 2*pi));
title('Small Change - Pre');
%---------------
subplot(2,3,3);
rh = rose2(mod(pre(catbin == 3)+2*pi, 2*pi));
title('Away - Pre');
%---------------
subplot(2,3,4);
rh = rose2(mod(post(catbin == 1)+2*pi, 2*pi));
title('Approach - Post Turn');
%---------------
subplot(2,3,5);
rh = rose2(mod(post(catbin == 2)+2*pi, 2*pi));
title('Small Change - Post Turn');
%---------------
subplot(2,3,6);
rh = rose2(mod(post(catbin == 3)+2*pi, 2*pi));
title('Away - Post Turn');

% Figure 3,4 - Turning Magnitude Scatterplot and Histograms
% -------------------------------------
plotblue = [.3 .6 1];
scat_fig = figure;
% Distance decreases - approaching
pre_angles = mod(sniffData.turnTrig_preTurnHeadings(catbin == 1)+(2*pi), 2*pi);
post_angles = mod(sniffData.turnTrig_postTurnHeadings(catbin == 1)+(2*pi), 2*pi);
plot(pre_angles, post_angles, '^g'); hold on;
diff_fig = figure; 
labelah = axes('Position', [ 0 0 1 1], 'Visible', 'off'); 
text(.3, .1, ' Turning Magnitudes', 'FontSize', 18);
subplot(2,3,1);
ang_diff = circ_dist(post_angles, pre_angles);
rose2(ang_diff);
title('Intersniff Approach');
%[~, ang_diff] = rotationDirection(pre_angles, post_angles);
med_diff = median(abs(ang_diff))
% Magnitude histograms
subplot(2,3,4); [N, cent] = hist(abs(ang_diff)); bar(cent, N, 'hist'); 
xlim([0 pi]); hold on; plot(gca, [med_diff med_diff], [0 max(N) + 10], 'k--');
text(med_diff, max(N)+3, ['median = ' num2str(med_diff)]);
set(gca, 'Xtick', [0 pi/2 pi], 'XTickLabel', {'0', 'p/2', 'p'}, 'FontName', 'symbol');
ha = findobj(gca,'Type','patch'); set(ha, 'FaceColor', plotblue);
xlabel('Turn Magnitude (radians)', 'FontName', 'Helvetica'); ylabel('Count', 'FontName', 'Helvetica');

% Distance spanning zero
pre_angles = mod(sniffData.turnTrig_preTurnHeadings(catbin == 2)+(2*pi), 2*pi);
post_angles = mod(sniffData.turnTrig_postTurnHeadings(catbin == 2)+(2*pi), 2*pi);
figure(scat_fig); plot(pre_angles, post_angles, 'xk');
figure(diff_fig); subplot(2,3,2);
ang_diff = circ_dist(post_angles, pre_angles);
rose2(ang_diff);
title('Intersniff Small Change');
%[~, ang_diff] = rotationDirection(pre_angles, post_angles);
med_diff = median(abs(ang_diff))
% Magnitude histograms
subplot(2,3,5); [N, cent] = hist(abs(ang_diff)); bar(cent, N, 'hist'); 
xlim([0 pi]); hold on; plot(gca, [med_diff med_diff], [0 max(N) + 10], 'k--');
text(med_diff, max(N)+3, ['median = ' num2str(med_diff)]);
set(gca, 'Xtick', [0 pi/2 pi], 'XTickLabel', {'0', 'p/2', 'p'}, 'FontName', 'symbol');
ha = findobj(gca,'Type','patch'); set(ha, 'FaceColor', plotblue);
xlabel('Turn Magnitude (radians)', 'FontName', 'Helvetica'); ylabel('Count', 'FontName', 'Helvetica');
title('Turns', 'FontName', 'Helvetica');
% Distance increases - away
pre_angles = mod(sniffData.turnTrig_preTurnHeadings(catbin == 3)+(2*pi), 2*pi);
post_angles = mod(sniffData.turnTrig_postTurnHeadings(catbin == 3)+(2*pi), 2*pi);
figure(scat_fig); plot(pre_angles, post_angles, 'or');
xlim([0 2*pi]); ylim([0 2*pi]);
figure(diff_fig); subplot(2,3,3);
ang_diff = circ_dist(post_angles, pre_angles);
rose2(ang_diff);
title('Intersniff Away');
%[~, ang_diff] = rotationDirection(pre_angles, post_angles);
med_diff = median(abs(ang_diff))
% Magnitude histograms
subplot(2,3,6); [N, cent] = hist(abs(ang_diff)); bar(cent, N, 'hist'); 
xlim([0 pi]); hold on; plot(gca, [med_diff med_diff], [0 max(N) + 10], 'k--');
text(med_diff, max(N)+3, ['median = ' num2str(med_diff)]);
set(gca, 'Xtick', [0 pi/2 pi], 'XTickLabel', {'0', 'p/2', 'p'}, 'FontName', 'symbol');
ha = findobj(gca,'Type','patch'); set(ha, 'FaceColor', plotblue);
xlabel('Turn Magnitude (radians)', 'FontName', 'Helvetica'); ylabel('Count', 'FontName', 'Helvetica');
%axes(labelah);
%text(.3, .1, 'Turning Magnitudes');

%% Figure 5 - Dependence of mean turning magnitude on the distance differences and the 
%             absolute position of the mouse.
corder = {[.9 .4 .4], [1 0 0], [0 0 0]};
edges_dd = -10:4:10;
edges_pos = -14:4:14;
post = convertMouseHeadingAngles(sniffData.turnTrig_postTurnHeadings, sniffData.turnTrig_postTurnPos, sniffData.turnTrig_postTurnDirToTrail);
pre = convertMouseHeadingAngles(sniffData.turnTrig_preTurnHeadings, sniffData.turnTrig_preTurnPos, sniffData.turnTrig_preTurnDirToTrail);
ang_diff = circ_dist(post, pre);
turn_mag = abs(ang_diff);
med_mags_dd = NaN*zeros(length(edges_dd),3);
mean_mags_dd = NaN*zeros(length(edges_dd),3);
ste_mags_dd = NaN*zeros(length(edges_dd),3);
med_mags_pos = NaN*zeros(length(edges_pos),3);
mean_mags_pos = NaN*zeros(length(edges_dd),3);
ste_mags_pos = NaN*zeros(length(edges_dd),3);
all_mags_dd = {};
for jj = 1:3
    [N, bin] = histc(sniffData.turnTrig_preTurnDistDiff(:,jj), edges_dd);
    for ii = 1:length(N)
        mags = turn_mag(bin==ii);
        med_mags_dd(ii, jj) = median(mags);
        mean_mags_dd(ii, jj) = mean(mags);
        ste_mags_dd(ii, jj) = std(mags)./sqrt(numel(mags));
        all_mags_dd{ii, jj} = mags;
    end
    
    [N, bin] = histc(sniffData.turnTrig_sniffPos(:,jj+1), edges_pos);
    for ii = 1:length(N)
        mags = turn_mag(bin==ii);
        med_mags_pos(ii, jj) = median(mags);
        mean_mags_pos(ii, jj) = mean(mags);
        ste_mags_pos(ii, jj) = std(mags)./sqrt(numel(mags));
        %all_mags_dd{ii, jj} = mags;
    end
end

% Two way ANOVA on the turning magnitudes
all_mags_dd = all_mags_dd(1:end-1, :);
buildANOVA2Matrix;
[p, table, stats] = anova2(testM, minelem);
[c,m,h] = multcompare(stats, 'ctype', 'bonferroni');


figure;
subplot(2,2,1);
plot(edges_dd', med_mags_dd, 'k-o');
xlabel('\Delta ND (mm)'); ylabel('Median Turn Mag'); xlim([-10 10]);
subplot(2,2,2); hold on;
x = edges_dd(1:end-1) + (edges_dd(2)-edges_dd(1))/2;
for ii = 1:3
    plot_err_poly(gca, x', mean_mags_dd(1:end-1,ii), ste_mags_dd(1:end-1,ii), corder{ii}, (corder{ii}+[1 1 1])/2);
end
xlabel('\Delta ND (mm)'); ylabel('Mean Turn Mag'); xlim([-10 10]);
subplot(2,2,3);
plot(edges_pos', med_mags_pos,'k-o');
xlabel('ND (mm)'); ylabel('Median Turn Mag'); xlim([-15 15]);
subplot(2,2,4); hold on;
x = edges_pos(1:end-1) + (edges_pos(2)-edges_pos(1))/2;
for ii = 1:3
    plot_err_poly(gca, x', mean_mags_pos(1:end-1,ii), ste_mags_pos(1:end-1,ii), corder{ii}, (corder{ii}+[1 1 1])/2);
end
xlabel('ND (mm)'); ylabel('Mean Turn Mag'); xlim([-15 15]);
%plot(sniffData.turnTrig_preTurnDistDiff(:,3), turn_mag, '.k', 'MarkerSize', 6);
%subplot(1,3,3);
%magM = cell2mat_2Dunequal(all_mags{:,3});
%boxplot(magM, 'notch', 'on');

%% Figure 6 - Plotting measures of turning relative to sniff location
% rather than intersniff differences

turnDirections = sniffData.turnDirections./abs(sniffData.turnDirections); %convert to signed 1's.
turnPos = sniffData.sniffPos(:,end);
turnsLeft = turnDirections.*turnPos;
edges = -15:2:15;
turnp = NaN*zeros(numel(edges), size(sniffData.sniffPos,2)-1);
for jj = 2:size(sniffData.sniffPos,2)
    [N,bin] = histc(sniffData.sniffPos(:,jj), edges);
    nbins = length(N);
    turn_counts = zeros(1,nbins);
    turnTowards= zeros(1,nbins);
    %turnTowardsProp = zeros(1,nbins);
    for ii = 1:nbins
        turn_counts(ii) = sum(sniffData.bturn(bin == ii));
        temp = turnDirections(bin == ii); % turns after each concentration change
        turnTowards(ii) = sum(temp(temp == 1)); %The number of turns towards the trail
        turnTowardsProp(ii, jj-1) = turnTowards(ii)/numel(temp(~isnan(temp)));
        temp = turnsLeft(bin == ii);
        nTurnsLeft = sum(temp > 0);
        pTurnsLeft(ii, jj-1) = nTurnsLeft./numel(temp(~isnan(temp)));
    end
    turnp(:,jj-1) = turn_counts(:)./N(:);
end

figure;
subplot(1, 3, 1);
stairs(edges',turnp);
xlabel('Distance from Trail (mm)');
ylabel('Turning Probability');
xlim([-15 15]);
ylim([0 .5]);

subplot(1,3,2);
stairs(edges', turnTowardsProp);
xlabel('Distance from Trail (mm)');
ylabel('Prop Turns Towards Trail');
legend({'Sniff -2', 'Sniff -1','Sniff 0'});
xlim([-15 15]);

subplot(1,3,3);
stairs(edges', pTurnsLeft);
xlabel('Distance from Trail (mm)');
ylabel('Prop Turns Leftwards');
legend({'Sniff -2', 'Sniff -1','Sniff 0'});
xlim([-15 15]);

%% Figure 7 - Turning Magnitudes as a function of sniff position

