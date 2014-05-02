% Analysis of the sniff phases of a trail crossing
% 
% Script, assumes that exp structure is loaded

mm_conv = 1.16; %mm/px linear
thresh_dist = 20;
traj_wind = -15:60;
upn = 4;
mov_thresh = 50; %mm/sec
% Empty variables, to start with
exp_sniffDists = []; exp_followDists = [];
exp_traj = []; exp_dir = []; exp_relSniffFrame = []; exp_relSniffDists=[];
exp_sniff_phase = []; exp_res_dists = [];
%Lets loop through the trials
for ii = 1:length(exp.resp)
    exp.vids(ii).makePathsSkel();
    % Need to select segments based on the velocity of the movement - don't want to analyze stopped mice
    noseVel = exp.vids(ii).noseVel * mm_conv * exp.vids(ii).frameRate; %get the nose/body velocities
    nosePos = interpM(mm_conv * exp.vids(ii).nosePos, upn);
    bodyVel = exp.vids(ii).bodyVel(:,1) * mm_conv * exp.vids(ii).frameRate;
    noseVel_filt = interp(gaussianFilter(noseVel, 3, 'conv'), upn); %smoother versions - vels tend to be noisy
    bodyVel_filt = interp(gaussianFilter(bodyVel, 3, 'conv'), upn); 
    sniffphase = assignSniffPhase(exp.resp(ii).sniffVect);
    % select the frames to analyze - when the animal is following the trail
    % selection criteria
    odist = double(mm_conv * exp.vids(ii).orthogonalDistFromTrail(1:exp.vids(ii).nFrames, 1));
    allDists = interp(odist, upn);
    followAll = abs(allDists) <= thresh_dist; %considered following trail
    moving = noseVel_filt >= mov_thresh;   movingInds = find(moving);
    % Get the trajectories of the nose relative to the trail
    [nose_traj, dirs, wind, crossingInds] = noseTrajectories(exp.vids(ii), traj_wind); %upsample 
    nose_traj = interpM(nose_traj',upn); nose_traj = nose_traj';
    [crossingFrames, cinds] = intersect(crossingInds, movingInds); %now the frames when moving and crossing
    nose_traj = mm_conv * nose_traj(cinds, :); dirs = dirs(cinds); % convert and select
    %adjust the frame timings by the offset found through interpolation of positions 
    exp_wind = interp(traj_wind, upn);
    interp_wind = find(exp_wind >= 0 & exp_wind <= 1);
    [res_dists, adj_i] = min(abs(nose_traj(:, interp_wind)'));
    exp_t = interp(exp.vids(ii).times, upn); 
    adj_t = exp.vids(ii).times(crossingFrames) - (upn-adj_i)/(upn*exp.vids(ii).frameRate);
    crossing_resp_i = indicesFromTimes(exp.resp(ii).vidTime, adj_t);
    sp = sniffphase(crossing_resp_i);
    %phase_inds = exp.camTrig(ii).frameInds(crossingFrames + double(exp.camTrig(ii).vidStartFrame) - 1);
    %sp = sniffphase(phase_inds);
    %res_dists = exp.vids(ii).orthogonalDistFromTrail(crossingFrames, 1);
    exp_sniff_phase = cat(1, exp_sniff_phase, sp); exp_res_dists =  cat(1, exp_res_dists, res_dists');
end

% Circular histogram!
figure;
lh = rose(exp_sniff_phase,20);
set(lh,'Color', 'k', 'LineWidth',2);
x = get(lh,'Xdata');
y = get(lh,'Ydata');
g=patch(x,y,'b');

p_val = circ_rtest(exp_sniff_phase);
disp(['The p-value for the Raleigh test for uniformity is ' num2str(p_val)]);

% Look at the residual distances
figure;
hist(exp_res_dists,30)