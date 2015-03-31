% Analysis of sniffing triggered average trajectories
%
% Script, assumes that exp structure is loaded

mm_conv = .862; %mm/px linear
thresh_dist = 20;
traj_wind = -30:30; 
exp_sniffDists = []; exp_followDists = [];
exp_traj = []; exp_dir = []; exp_relSniffTime = []; exp_relSniffDists=[];
exp_nosePos = []; %nose position throughout the exp
ns = 8; %number of sniffs to analyze
upn = 10; %how much to upsample the nose position data, x per sample
long_wind = (min(traj_wind)*upn):(max(traj_wind)*upn);
cm = jet(ns); %color to plot the sniff scatter plot in
for ii = 1:length(exp.resp)
    exp.vids(ii).makePathsSkel();
    % Velocities
    noseVel = exp.vids(ii).noseVel * mm_conv * exp.vids(ii).frameRate; %get the nose/body velocities
    nosePos = interpM(mm_conv * exp.vids(ii).nosePos, upn);
    bodyVel = exp.vids(ii).bodyVel(:,1) * mm_conv * exp.vids(ii).frameRate;
    noseVel_filt = interp(gaussianFilter(noseVel, 1.5, 'conv'), upn); %smoother versions - vels tend to look messy
    bodyVel_filt = interp(gaussianFilter(bodyVel, 1.5, 'conv'), upn);
    sniffphase = assignSniffPhase(exp.resp(ii).sniffVect); %this line doesn't work
  
    odist = double(mm_conv * exp.vids(ii).orthogonalDistFromTrail(1:exp.vids(ii).nFrames, 1));
    fdist = odist(abs(odist) < thresh_dist);
    exp_nosePos = cat(1, exp_nosePos, fdist(:));
    allDists = interp(odist, upn);
    % Get the trajectories of the nose relative to the trail
    [turningInds, dirs, nose_traj] = exp.vids(ii).findFollowingTurns([], 1, thresh_dist, traj_wind);
    nose_traj = interpMatrix(nose_traj,upn); nose_traj = nose_traj'; %upsample - interpolate the position
    %[crossingFrames, cinds] = intersect(turningInds, movingInds); %now the frames when moving and crossing
    %nose_traj = mm_conv * nose_traj(cinds, :); dirs = dirs(cinds);
    exp_traj = cat(1, exp_traj, nose_traj); exp_dir = cat(1,exp_dir, dirs);
    %adjust the frame timings by the offset found through interpolation of positions 
    exp_wind = interp(traj_wind, upn);
    interp_wind = find(exp_wind >= 0 & exp_wind <= 1);
    %[res_dists, adj_i] = min(abs(nose_traj(:, interp_wind)'));
    exp_t = interp(exp.vids(ii).times, upn); % longer time vector
    
    
    %crossing_resp_i = indicesFromTimes(exp.resp(ii).vidTime, exp_t); % what?
    %sp = sniffphase(crossing_resp_i);
    sniffTimes = exp.resp(ii).vidSniffTimes(exp.resp(ii).vidSniffs);
    newTurningInds = (turningInds-1)*upn + 1;
    newTurningTimes = exp_t(newTurningInds);
    %newCrossingInds = (crossingFrames-1)*upn + 1 - (upn-adj_i');
    % want to build the distribution of sniff times and positions relative to crossing
    relSniffTime = NaN*zeros(length(newTurningInds), ns);
    relSniffDists = NaN*zeros(length(newTurningInds), ns);
    for jj=1:length(newTurningInds) 
        sniffBefore = find(sniffTimes <= newTurningTimes(jj),4, 'last');
        if ~isempty(sniffBefore)
            si = sniffBefore(1):(sniffBefore(1)+ns-1);
            if sum(si > length(sniffTimes)) %check for array overrun 
                break;
            end
            relSniffTime(jj,:) = sniffTimes(si) - newTurningTimes(jj);
            relSniffDists(jj,:) = allDists(indicesFromTimes(exp_t,sniffTimes(si)));
        end
        
    end
    exp_relSniffTime = cat(1, exp_relSniffTime, relSniffTime*1000);
    exp_relSniffDists = cat(1, exp_relSniffDists, relSniffDists);
    
end
exp_wind = exp_wind / exp.vids(1).frameRate * 1000; %convert to ms

%%
% We should do some calculation for the Khan et al (Bhalla), turn-triggered difference of sniff position
mean_pos = nanmean(exp_nosePos);
absDists = abs(exp_relSniffDists-mean_pos);
distChange = diff(absDists, 1, 2);
mean_distChange = nanmean(distChange);
label = {'-4to-3', '-3to-2', '-2to-1', '-1to1', '1to2','2to3', '3to4'};
dirs = [-1,1]; 
figure; 
for ii=1:length(dirs)
    curr_dir = exp_dir == dirs(ii);
    mean_distChange(ii,:) = nanmean(distChange(curr_dir,:));
end
bar(mean_distChange'); 
set(gca, 'XTickLabel', label)
xlabel('Sniff relative to turn');
ylabel('Distance Change from center of trail');
legend({'Rightward Turn', 'Leftward Turn'});

%% One complicated figure.
% Main axis - individual nose trajectories (distances from trail) in grey, averaged in black.  Sniff
% positions and times as colored points overlaid to visualize when and where the animals sniffs.
% Top axis - the marginal distribution of sniffing times
% Right axis - the marginal distribution of the sniffing positions.  
f1= figure;
main_ax = axes('Position', [.1 .1 .6 .6]); hold on;
right_ax = axes('Position', [.8 .1 .15 .6]); hold on;
top_ax = axes('Position', [.1 .8 .6 .15]); hold on;
dirs = [-1,1]; dirc = {[.4 0 .4], [0 .4 .4]};
return_dir = [-1, 1];
hist_tbin = long_wind;
sel_time = [0 50];
for ii=1:length(dirs)
    curr_dir = exp_dir == dirs(ii);
    mean_traj = nanmean(exp_traj(curr_dir,:));
    mean_sniff_dist = nanmean(exp_relSniffDists(curr_dir, :));
    mean_sniff_time = nanmean(exp_relSniffTime(curr_dir, :));
    % select trials to see that have a direction change within 500ms
    sel_traj = exp_traj(curr_dir,:); %subset in the correct direction
    
    sel = false(size(sel_traj,1),1); %sel(1:40:end) = true;
    %hold on; plot(main_ax, exp_wind, sel_traj(sel,:), 'Color', dirc{ii}, 'Linewidth',.5);
    hold on; plot(main_ax, exp_wind, mean_traj, 'k', 'Linewidth',2);
    sel_relSniffTime = exp_relSniffTime(:,ii);
    % plot the sniff positions in a different color relative to the direction change
    sel = false(size(curr_dir,1),1); %sel(1:10:end) = true;
    sel = curr_dir & sel;
    for jj=1:ns 
        plot(main_ax, exp_relSniffTime(sel,jj), exp_relSniffDists(sel,jj), '.', 'Color', cm(jj,:));
        plot(main_ax, mean_sniff_time(jj), mean_sniff_dist(jj), '.' , 'Color', cm(jj,:), 'MarkerSize', 12);
    end
end
plot(main_ax, exp_wind, zeros(size(exp_wind)), '--k', 'LineWidth', 1);
% for ii=1:ns
%    plot(main_ax, exp_relSniffTime(:,ii), exp_relSniffDists(:,ii), '.', 'Color', cm(ii,:));
% end
hist_bin = exp_wind(1:2:length(exp_wind));
topy = histc(exp_relSniffTime(:), hist_bin);
bar(top_ax, hist_bin, topy, 'r', 'EdgeColor', 'r', 'BarWidth', 1);
%topy = histc(exp_relSniffTime(:), hist_tbin);
%bar(top_ax, hist_tbin, topy, 'r');

dist_lim = [-20 20];
[sd, x] = hist(exp_relSniffDists(:), 100); 
barh(right_ax, x,sd, 'r', 'BarWidth', 1);
set(main_ax, 'Xlim', [exp_wind(1) exp_wind(end)], 'ylim',dist_lim,'YDir', 'reverse');
set(top_ax, 'xlim', [exp_wind(1) exp_wind(end)]);
set(right_ax, 'ylim', dist_lim, 'Ydir', 'reverse');
xlabel(main_ax, 'Time relative to crossing');
ylabel(main_ax, 'Distance from trail center (mm)');

%% Also want to test the uniformity of sniffing times using a Chi-squared test
sniffi = exp_relSniffTime(:) >= -200 & exp_relSniffTime(:) <= 200;
edges = linspace(-200,200, 40);
centers = edges(1:end-1)+ 5;
obs = histc(exp_relSniffTime(sniffi),edges); 
gen = -200 + 400 .* rand(2000,1);
geny = histc(gen, edges);
expected = sum(obs)./(length(edges)-1) * ones(length(edges)-1, 1);
expected_cdf = cumsum(expected);
expected_rand = sum(geny)./(length(edges)-1) * ones(length(edges)-1, 1);
expected_rand_cdf = cumsum(expected_rand);
[h, p, stats] = chi2gof(obs(1:end-1), 'expected', expected, 'nparams', 0);  
disp(['Chi square test for uniformity: h = ' num2str(h) '   p = ' num2str(p)]);
figure; bar(edges, obs, 'histc');
hold on; plot([-200 200], [expected(1) expected(end)], '--k');
[h, p, stats] = chi2gof(geny(1:end-1), 'expected', expected_rand, 'nparams', 0);  
disp(['Chi square test for uniformity (random data): h = ' num2str(h) '   p = ' num2str(p)]);
figure; bar(edges, geny, 'histc');
hold on; plot([-200 200], [expected_rand(1) expected_rand(end)], '--k');
% Let's do the analogous KS test on it.
expected_cdf = cumsum(expected); expected_cdf = [expected_cdf; expected_cdf(end)]; %create an expected CDF
expected_rand_cdf = cumsum(expected_rand); expected_rand_cdf = [expected_rand_cdf; expected_rand_cdf(end)]; %create an expected CDF
obs = obs(1:end-1); obs_cdf = cumsum(obs(:));
geny = geny(1:end-1); geny_cdf = cumsum(geny(:));
[h,p,ksstat,cv] = kstest(exp_relSniffTime(sniffi), [edges(:), expected_cdf]);
disp(['KS test for uniformity (random data): h = ' num2str(h) '   p = ' num2str(p)]);
[h2,p2,ksstat2,cv2] = kstest(gen, [edges(:), expected_rand_cdf]);
disp(['KS test for uniformity (random data): h = ' num2str(h2) '   p = ' num2str(p2)]);


%% Trying to recreate the Khan et al result that before a crossing/turning you are different distances
%  from the path.  

