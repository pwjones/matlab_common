% Analysis of sniffing triggered average trajectories
%
% Script, assumes that exp structure is loaded

mm_conv = 1.16; %mm/px linear
thresh_dist = 20;
traj_wind = -15:60;
exp_sniffDists = []; exp_followDists = [];
exp_traj = []; exp_dir = []; exp_relSniffFrame = []; exp_relSniffDists=[];
ns = 8; %number of sniffs to analyze
cm = jet(ns); %color to plot the sniff scatter plot in
for ii = 1:length(exp.resp)
    exp.vids(ii).makePathsSkel();
    % Velocities
    noseVel = exp.vids(ii).noseVel * mm_conv * exp.vids(ii).frameRate; %get the nose/body velocities
    nosePos = mm_conv * exp.vids(ii).nosePos ;
    bodyVel = exp.vids(ii).bodyVel(:,1) * mm_conv * exp.vids(ii).frameRate;
    noseVel_filt = gaussianFilter(noseVel, 3, 'conv'); %smoother versions - vels tend to look messy
    bodyVel_filt = gaussianFilter(bodyVel, 3, 'conv'); 
    % select the frames to analyze - when the animal is following the trail
    % selection criteria
    allDists = mm_conv * exp.vids(ii).orthogonalDistFromTrail(1:exp.vids(ii).nFrames, 1);
    followAll = abs(allDists) <= thresh_dist;
    moving = noseVel_filt >= 50;
    movingInds = find(moving);
    [nose_traj, dirs, wind, crossingInds] = noseTrajectories(exp.vids(ii), traj_wind);
    %nose_traj = interpM(nose_traj',4); nose_traj = nose_traj';
    [crossingFrames, cinds] = intersect(crossingInds, movingInds); %now the frames when moving and crossing
    nose_traj = mm_conv * nose_traj(cinds, :); dirs = dirs(cinds);
    sniffFrames = exp.resp(ii).sniffFrames(exp.resp(ii).vidSniffs);
    exp_traj = cat(1, exp_traj, nose_traj); exp_dir = cat(1,exp_dir, dirs);
    % want to build the distribution of sniff times and positions relative to crossing
    relSniffFrame = NaN*zeros(length(crossingFrames), ns);
    relSniffDists = NaN*zeros(length(crossingFrames), ns);
    for jj=1:length(crossingFrames) 
        sniffBefore = find(sniffFrames <= crossingFrames(jj),4, 'last');
        if ~isempty(sniffBefore)
            si = sniffBefore(1):(sniffBefore(1)+ns-1);
            if sum(si > length(sniffFrames)) %check for array overrun 
                break;
            end
            relSniffFrame(jj,:) = sniffFrames(si) - crossingFrames(jj);
            relSniffDists(jj,:) = allDists(sniffFrames(si));
        end
    end
    exp_relSniffFrame = cat(1, exp_relSniffFrame, relSniffFrame);
    exp_relSniffDists = cat(1, exp_relSniffDists, relSniffDists);
    
end


f1= figure;
main_ax = axes('Position', [.1 .1 .6 .6]); hold on;
right_ax = axes('Position', [.8 .1 .15 .6]); hold on;
top_ax = axes('Position', [.1 .8 .6 .15]); hold on;
dirs = [0,1];
return_dir = [-1, 1];
traj_wind = traj_wind / exp.vids(1).frameRate * 1000; %convert to ms
sel_time = [0 500];
for ii=1:length(dirs)
    curr_dir = exp_dir == dirs(ii);
    mean_traj = nanmean(exp_traj(curr_dir,:));
    % select trials to see that have a direction change within 500ms
    sel_traj = exp_traj(curr_dir,:);
    sel_ti = traj_wind >= sel_time(1) & traj_wind <= sel_time(2);
    traj_prime = diff(sel_traj,1,2);
    with_return = false(size(sel_traj,1),1);
    for jj=1:size(sel_traj,1)
        test = traj_prime(jj,sel_ti) .* return_dir(ii);
        if sum(test>0) %there are any trajectories of the return dir
            with_return(jj) = 1;
        end
    end
    sel = false(size(sel_traj,1),1); sel(1:8:end) = 1;
    hold on; plot(main_ax, traj_wind, sel_traj(with_return & sel,:), 'Color', [.4 .4 .4], 'Linewidth',.5);
    hold on; plot(main_ax, traj_wind, mean_traj, 'k', 'Linewidth',2);
end
plot(main_ax, traj_wind, zeros(size(traj_wind)), '--k', 'LineWidth', .5);
for ii=1:ns
    plot(main_ax, exp_relSniffFrame(:,ii) ./ exp.vids(1).frameRate * 1000, exp_relSniffDists(:,ii), '.', 'Color', cm(ii,:));
end
topy = histc(exp_relSniffFrame(:) ./exp.vids(1).frameRate * 1000, traj_wind);
bar(top_ax, traj_wind, topy, 'r');
[sd, x] = hist(exp_relSniffDists(:), 40); 
barh(right_ax, x,sd, 'r');
set(main_ax, 'Xlim', [traj_wind(1) 1000], 'ylim', [-40 40]);
set(top_ax, 'xlim', [traj_wind(1) 1000]);
set(right_ax, 'ylim', [-40 40]);
xlabel(main_ax, 'Time relative to crossing');
ylabel(main_ax, 'Distance from trail center (mm)');