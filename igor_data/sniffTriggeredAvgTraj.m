% Analysis of sniffing triggered average trajectories
%
% Script, assumes that exp structure is loaded

mm_conv = 1.16; %mm/px linear
thresh_dist = 20;
traj_wind = -30:60;
exp_sniffDists = []; exp_followDists = [];
exp_traj = []; exp_dir = []; exp_relSniffTime = []; exp_relSniffDists=[];
ns = 8; %number of sniffs to analyze
upn = 4; %how much to upsample the nose position data, x per sample
long_wind = (min(traj_wind)*upn):(max(traj_wind)*upn);
cm = jet(ns); %color to plot the sniff scatter plot in
for ii = 1:length(exp.resp)
    exp.vids(ii).makePathsSkel();
    % Velocities
    noseVel = exp.vids(ii).noseVel * mm_conv * exp.vids(ii).frameRate; %get the nose/body velocities
    nosePos = interpM(mm_conv * exp.vids(ii).nosePos, upn);
    bodyVel = exp.vids(ii).bodyVel(:,1) * mm_conv * exp.vids(ii).frameRate;
    noseVel_filt = interp(gaussianFilter(noseVel, 3, 'conv'), upn); %smoother versions - vels tend to look messy
    bodyVel_filt = interp(gaussianFilter(bodyVel, 3, 'conv'), upn);
    sniffphase = assignSniffPhase(exp.resp(ii).sniffVect);
    % select the frames to analyze - when the animal is following the trail
    % selection criteria
    odist = double(mm_conv * exp.vids(ii).orthogonalDistFromTrail(1:exp.vids(ii).nFrames, 1));
    allDists = interp(odist, upn);
    followAll = abs(allDists) <= thresh_dist;
    moving = noseVel_filt >= 50;    movingInds = find(moving);
    % Get the trajectories of the nose relative to the trail
    [nose_traj, dirs, wind, crossingInds] = noseTrajectories(exp.vids(ii), traj_wind); %upsample 
    nose_traj = interpM(nose_traj',upn); nose_traj = nose_traj';
    [crossingFrames, cinds] = intersect(crossingInds, movingInds); %now the frames when moving and crossing
    nose_traj = mm_conv * nose_traj(cinds, :); dirs = dirs(cinds);
    exp_traj = cat(1, exp_traj, nose_traj); exp_dir = cat(1,exp_dir, dirs);
    %adjust the frame timings by the offset found through interpolation of positions 
    exp_wind = interp(traj_wind, upn);
    interp_wind = find(exp_wind >= 0 & exp_wind <= 1);
    [res_dists, adj_i] = min(abs(nose_traj(:, interp_wind)'));
    exp_t = interp(exp.vids(ii).times, upn); % longer time vector
    adj_t = exp.vids(ii).times(crossingFrames) - (upn-adj_i)/(upn*exp.vids(ii).frameRate); %adjusted times
    crossing_resp_i = indicesFromTimes(exp.resp(ii).vidTime, adj_t);
    sp = sniffphase(crossing_resp_i);
    sniffTimes = exp.resp(ii).vidSniffTimes(exp.resp(ii).vidSniffs);
    newCrossingInds = (crossingFrames-1)*upn + 1 - (upn-adj_i');
    % want to build the distribution of sniff times and positions relative to crossing
    relSniffTime = NaN*zeros(length(crossingFrames), ns);
    relSniffDists = NaN*zeros(length(crossingFrames), ns);
    for jj=1:length(crossingFrames) 
        sniffBefore = find(sniffTimes <= adj_t(jj),4, 'last');
        if ~isempty(sniffBefore)
            si = sniffBefore(1):(sniffBefore(1)+ns-1);
            if sum(si > length(sniffFrames)) %check for array overrun 
                break;
            end
            relSniffTime(jj,:) = sniffTimes(si) - adj_t(jj);
            relSniffDists(jj,:) = allDists(indicesFromTimes(exp_t,sniffTimes(si)));
        end
        
    end
    exp_relSniffTime = cat(1, exp_relSniffTime, relSniffTime*1000);
    exp_relSniffDists = cat(1, exp_relSniffDists, relSniffDists);
    
end

% One complicated figure.
% Main axis - individual nose trajectories (distances from trail) in grey, averaged in black.  Sniff
% positions and times as colored points overlaid to visualize when and where the animals sniffs.
% Top axis - the marginal distribution of sniffing times
% Right axis - the marginal distribution of the sniffing positions.  
f1= figure;
main_ax = axes('Position', [.1 .1 .6 .6]); hold on;
right_ax = axes('Position', [.8 .1 .15 .6]); hold on;
top_ax = axes('Position', [.1 .8 .6 .15]); hold on;
dirs = [0,1];
return_dir = [-1, 1];
exp_wind = exp_wind / exp.vids(1).frameRate * 1000; %convert to ms
sel_time = [0 500];
for ii=1:length(dirs)
    curr_dir = exp_dir == dirs(ii);
    mean_traj = nanmean(exp_traj(curr_dir,:));
    % select trials to see that have a direction change within 500ms
    sel_traj = exp_traj(curr_dir,:);
    sel_ti = exp_wind >= sel_time(1) & exp_wind <= sel_time(2);
    traj_prime = diff(sel_traj,1,2);
    with_return = false(size(sel_traj,1),1);
    for jj=1:size(sel_traj,1)
        test = traj_prime(jj,sel_ti) .* return_dir(ii);
        if sum(test>0) %there are any trajectories of the return dir
            with_return(jj) = 1;
        end
    end
    sel = false(size(sel_traj,1),1); sel(1:8:end) = 1;
    hold on; plot(main_ax, exp_wind, sel_traj(with_return & sel,:), 'Color', [.4 .4 .4], 'Linewidth',.5);
    hold on; plot(main_ax, exp_wind, mean_traj, 'k', 'Linewidth',2);
end
plot(main_ax, exp_wind, zeros(size(exp_wind)), '--k', 'LineWidth', .5);
for ii=1:ns
    plot(main_ax, exp_relSniffTime(:,ii), exp_relSniffDists(:,ii), '.', 'Color', cm(ii,:));
end
topy = histc(exp_relSniffTime(:), exp_wind);
bar(top_ax, exp_wind, topy, 'r');
[sd, x] = hist(exp_relSniffDists(:), 40); 
barh(right_ax, x,sd, 'r');
set(main_ax, 'Xlim', [exp_wind(1) 1000], 'ylim', [-40 40]);
set(top_ax, 'xlim', [exp_wind(1) 1000]);
set(right_ax, 'ylim', [-40 40]);
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
