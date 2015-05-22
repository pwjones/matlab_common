% Analysis of sniffing triggered average trajectories
% Parsing the trajectories based on the intersniff distance change.
% Script, assumes that exp structure is loaded
dbg=0;
mm_conv = .862; %mm/px linear
thresh_dist = 20;
traj_wind = -5:15; 
exp_sniffDists = []; exp_followDists = [];
exp_traj = []; exp_dir = []; exp_relSniffTime = []; exp_relSniffDists=[];
exp_nosePos = []; %nose position throughout the exp
upn = 10; %how much to upsample the nose position data, x per sample
turn_wind = .08; %s, time after a sniff to look for turns
long_wind = (min(traj_wind)*upn):(max(traj_wind)*upn);

% data structures for the turning probabilities
dist_diffs = []; bturn = []; turnDirections = []; mouseHeading = [];
allPreTurnHeadings = []; allTurnTrajectories = []; allTurnDirs = []; 
allPostTurnHeadings = []; allPreTurnDistDiff = []; turnLag = [];
allSniffHeadings = [];
for ii = 1:length(exp.resp)
    
    exp.vids(ii).makePathsSkel();
    % Velocities
    noseVel = exp.vids(ii).noseVel * mm_conv * exp.vids(ii).frameRate; %get the nose/body velocities
    nosePos = interpM(mm_conv * exp.vids(ii).nosePos, upn);
    bodyVel = exp.vids(ii).bodyVel(:,1) * mm_conv * exp.vids(ii).frameRate;
    noseVel_filt = interp(gaussianFilter(noseVel, 1.5, 'conv'), upn); %smoother versions - vels tend to look messy
    bodyVel_filt = interp(gaussianFilter(bodyVel, 1.5, 'conv'), upn);
    %sniffphase = assignSniffPhase(exp.resp(ii).sniffVect); %this line doesn't work
    dt = exp.vids(ii).times(2) - exp.vids(ii).times(1); 
    % Let's gather the sniff information
    sniffFrames = exp.resp(ii).sniffFrames(exp.resp(ii).vidSniffs);
    sniffFrames_upsamp = (sniffFrames -1)*upn + 1;
    
    %The nose trajectories and following segments
    [odist, heading]= exp.vids(ii).orthogonalDistFromTrail(1:exp.vids(ii).nFrames, 1);
    odist = double(mm_conv * odist); 
    odist_upsamp = interp(odist, upn);
    fsegs = exp.vids(ii).getFollowingSegments([], 1, thresh_dist);
    
    % Get the turning information - these are turns that happen within the
    % thresh_dist from the trail, so it's already selected for following
    [turningInds, dirs, nose_traj] = exp.vids(ii).findFollowingTurns([], 1, thresh_dist, traj_wind);
    % Need to just exclude turns that are too early or late for proximal
    % events to be out of range
    sel = (turningInds + traj_wind(1)) > 0  & (turningInds - traj_wind(1) <= exp.vids(ii).nFrames);
    turningInds = turningInds(sel); dirs = dirs(sel); nose_traj = nose_traj(:,sel);
    
    % Finding the distance change between previous 2 sniffs
    dd = NaN*zeros(length(turningInds),1); 
    lag = NaN*zeros(length(turningInds),1);
    sniffi = NaN*ones(length(turningInds),1);
    for jj=1:length(turningInds) 
       preTS = find(sniffFrames < turningInds(jj), 2, 'last');
       if length(preTS) == 2 
           dd(jj) = diff(abs(odist(sniffFrames(preTS))));
           lag(jj) = turningInds(jj) - sniffFrames(preTS(2));
           sniffi(jj) = sniffFrames(preTS(2));
       end
    end
    sel = ~isnan(sniffi);
    turningInds = turningInds(sel); dirs = dirs(sel); nose_traj = nose_traj(:,sel);
    dd = dd(sel); lag = lag(sel); sniffi = sniffi(sel);
    
    tpos = exp.vids(ii).orthogonalDistFromTrail(turningInds+traj_wind(1),1);
    thead = exp.vids(ii).headingFromMotion_TrailRelative(turningInds+traj_wind(1), 1); %heading before turn
    theadPost = exp.vids(ii).headingFromMotion_TrailRelative(turningInds-traj_wind(1), 1); %heading after turn
    sniffhead = exp.vids(ii).headingFromMotion_TrailRelative(sniffi, 1); %heading at the sniff preceding
    
    tpos = double(mm_conv.*tpos);
    allTurnTrajectories = [allTurnTrajectories, nose_traj];
    allPreTurnHeadings = [allPreTurnHeadings; thead];
    allSniffHeadings = [allSniffHeadings; sniffhead];
    allPostTurnHeadings = [allPostTurnHeadings; theadPost];
    allTurnDirs = [allTurnDirs; dirs(:).*tpos(:)];
    allPreTurnDistDiff = [allPreTurnDistDiff; dd];
    

    % Its best if we do this following segment by segment, because if
    % we anneal them together and try to look at concentration differences
    % between sniffs then we'll get unrelated sniff concentration
    % differences.  
    if dbg
        figure; hold on;
    end
    followi = []; followi_upsamp = []; 
    for jj = 1:size(fsegs,1)
        range = fsegs(jj,1):fsegs(jj,2); %range of frames
        [turns, ~,ti] = intersect(range, turningInds); % selecting within range items
        sniffs = intersect(range, sniffFrames);
        turnDirs = dirs(ti); %turning directions
        
        for kk=1:(length(sniffs)-1)
            dist_diffs = [dist_diffs; abs(odist(sniffs(kk+1))) - abs(odist(sniffs(kk)))];
            mouseHeading = [mouseHeading; heading(sniffs(kk+1))]; %the mouse angle relative to the vector to trail
            frame_wind = sniffs(kk+1):(sniffs(kk+1) + floor(turn_wind/dt));
            [~,turni] = intersect(turns,frame_wind);
            if ~isempty(turni)
                didturn = 1;
                td = turnDirs(turni(1));
            else
                didturn = 0;
                td = NaN;
            end
            bturn = [bturn; didturn];
            % Defining the turn direction is towards/away from trail.
            % Previous turnDir is signed -1 for a rightward turn, 1 for a
            % leftward one.  Since the left side of the trail is negative
            % space, the positive product is toward the trail and neg is
            % away from the trail.
            turnDirections = [turnDirections; td*odist(sniffs(kk+1))];
        end
        
        followi = [followi, range];
        range_upsamp = ((fsegs(jj,1)-1)*upn+1):((fsegs(jj,2)-1)*upn+1);
        followi_upsamp = [followi_upsamp, range_upsamp];
        if dbg
            plot(exp.vids(ii).times(range), odist(range), '-b'); hold on;
        end
    end
    fdist = odist(followi); %following distances
    fdist_upsamp = odist_upsamp(followi_upsamp); %following distances
    exp_nosePos = cat(1, exp_nosePos, fdist(:)); 
    
    turningInds_upsamp = ((turningInds-1)*upn)+1;
    nose_traj = interpMatrix(nose_traj,upn); %upsample - interpolate the position
    nose_traj = nose_traj'; 
    exp_traj = cat(1, exp_traj, nose_traj); exp_dir = cat(1,exp_dir, dirs);
    exp_t = exp.vids(ii).times;
    t_upsamp = interp(exp.vids(ii).times, upn); % longer time vector
    

    if (dbg)
        %figure;
        %plot(t_upsamp(followi_upsamp), odist_upsamp(followi_upsamp));
        hold on;
        plot(exp_t(sniffAndFollow), odist(sniffAndFollow), 'o');
        plot(exp_t(turningInds), odist(turningInds),'x');
    end
    
end
%%
allPreTurnHeadings = mod(allPreTurnHeadings+2*pi, 2*pi); %make them all positive.
allPostTurnHeadings = mod(allPostTurnHeadings+2*pi, 2*pi);
turnDirections = turnDirections./abs(turnDirections); %convert to signed 1's.
edges = -9.5:9.5;
[N,bin] = histc(dist_diffs, edges);
nbins = length(N);
turn_counts = zeros(1,nbins); 
turnTowards= zeros(1,nbins);  
turnTowardsProp = zeros(1,nbins);
for ii = 1:nbins
    turn_counts(ii) = sum(bturn(bin == ii)); 
    temp = turnDirections(bin == ii); % turns after each concentration change
    turnTowards(ii) = sum(temp(temp == 1)); %The number of turns towards the trail
    turnTowardsProp(ii) = turnTowards(ii)/turn_counts(ii);
    %headings = mouseHeading(bin == ii);
    %meanHeadings(ii) = mean(headings);
    %stdHeadings(ii) = std(headings);
end

turnp = turn_counts(:)./N(:);
turnp(isnan(turnp)) = 0; %replace possible NaNs
 
figure;
subplot(2,2,1);
%bar(edges(1:end-1), N, 'histc');
bar(edges, N, 'histc');
xlabel('\Delta Distance from Trail (mm)');
ylabel('Counts');
xlim([-10 10]);

subplot(2,2,2);
%bar(edges(1:end-1),turnp, 'histc');
bar(edges,turnp, 'histc');
xlabel('\Delta Distance from Trail (mm)');
ylabel('Turning Probability');
xlim([-10 10]);
ylim([0 .5]);

subplot(2,2,3);
%bar(edges(1:end-1),turnTowardsProp, 'histc');
bar(edges, turnTowardsProp, 'histc');
xlabel('\Delta Distance from Trail (mm)');
ylabel('Prop Turns Towards Trail');
xlim([-10 10]);

%% Plotting Rose plots of before and after headings
c = [.3, .6, 1];
figure;
subplot(2,3,1);
cats = [-10 -2.5 2.5 10];
[Ncat, catbin] = histc(allPreTurnDistDiff, cats);
rh = rose2(mod(allSniffHeadings(catbin == 1)+2*pi, 2*pi));
title('Approach - Pre');
%---------------
subplot(2,3,2);
rh = rose2(mod(allSniffHeadings(catbin == 2)+2*pi, 2*pi));
title('Small Change - Pre');
%---------------
subplot(2,3,3);
rh = rose2(mod(allSniffHeadings(catbin == 3)+2*pi, 2*pi));
title('Away - Pre');
%---------------
subplot(2,3,4);
rh = rose2(mod(allPostTurnHeadings(catbin == 1)+2*pi, 2*pi));
title('Approach - Post Turn');
%---------------
subplot(2,3,5);
rh = rose2(mod(allPostTurnHeadings(catbin == 2)+2*pi, 2*pi));
title('Small Change - Post Turn');
%---------------
subplot(2,3,6);
rh = rose2(mod(allPostTurnHeadings(catbin == 3)+2*pi, 2*pi));
title('Away - Post Turn');
%---------------

%% Figures doing a scatter plot of the same idea - before vs after. Also, rose plots showing the magnitudes of heading differences
plotblue = [.3 .6 1];
scat_fig = figure;
% Distance decreases - approaching
pre_angles = mod(allPreTurnHeadings(catbin == 1)+(2*pi), 2*pi);
post_angles = mod(allPostTurnHeadings(catbin == 1)+(2*pi), 2*pi);
plot(pre_angles, post_angles, '^g'); hold on;
diff_fig = figure; 
labelah = axes('Position', [ 0 0 1 1], 'Visible', 'off'); 
text(labelah, .3, .1, ' Turning Magnitudes', 'FontSize', 18);
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
pre_angles = mod(allPreTurnHeadings(catbin == 2)+(2*pi), 2*pi);
post_angles = mod(allPostTurnHeadings(catbin == 2)+(2*pi), 2*pi);
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
pre_angles = mod(allPreTurnHeadings(catbin == 3)+(2*pi), 2*pi);
post_angles = mod(allPostTurnHeadings(catbin == 3)+(2*pi), 2*pi);
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
text(labelah, .3, .1, 'Turning Magnitudes');