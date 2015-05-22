% Analysis of sniffing triggered average trajectories
% Parsing the trajectories based on the intersniff distance change.
% Script, assumes that exp structure is loaded
dbg=0;
mm_conv = .862; %mm/px linear
thresh_dist = 20;
traj_wind = -30:30; 
exp_sniffDists = []; exp_followDists = [];
exp_traj = []; exp_dir = []; exp_relSniffTime = []; exp_relSniffDists=[];
exp_nosePos = []; %nose position throughout the exp
ns = 8; %number of sniffs to analyze
upn = 10; %how much to upsample the nose position data, x per sample
turn_wind = .080; %s, time after a sniff to look for turns
long_wind = (min(traj_wind)*upn):(max(traj_wind)*upn);
cm = jet(ns); %color to plot the sniff scatter plot in
% data structures for the turning probabilities
dist_diffs = []; bturn = [];
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
    odist = double(mm_conv * exp.vids(ii).orthogonalDistFromTrail(1:exp.vids(ii).nFrames, 1));
    odist_upsamp = interp(odist, upn);
    fsegs = exp.vids(ii).getFollowingSegments([], 1, thresh_dist);
    
    % Get the turning information 
    [turningInds, dirs, nose_traj] = exp.vids(ii).findFollowingTurns([], 1, thresh_dist, traj_wind);
    
    % So, its best if we do this following segment by segment, because if
    % we anneal them together and try to look at concentration differences
    % between sniffs then we'll get unrelated sniff concentration
    % differences.  
    if dbg
        figure; hold on;
    end
    followi = []; followi_upsamp = [];
    for jj = 1:size(fsegs,1)
        range = fsegs(jj,1):fsegs(jj,2); %range of frames
        turns = intersect(range, turningInds);
        sniffs = intersect(range, sniffFrames);
        
        for kk=1:(length(sniffs)-1)
            dist_diffs = [dist_diffs; abs(odist(sniffs(kk+1))) - abs(odist(sniffs(kk)))];
            frame_wind = sniffs(kk+1):(sniffs(kk+1) + floor(turn_wind/dt));
            turni = intersect(turns,frame_wind);
            if ~isempty(turni)
                didturn = 1;
            else
                didturn = 0;
            end
            bturn = [bturn; didturn];
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
edges = -9.5:1:9.5;
[N,edges,bin] = histcounts(dist_diffs, edges);
nbins = length(N);
turn_counts = zeros(1,nbins);
for ii = 1:nbins
    turn_counts(ii) = sum(bturn(bin == ii));
end
turnp = turn_counts./N;
turnp(isnan(turnp)) = 0; %replace possible NaNs
figure; 
bar(edges(1:end-1),turnp, 'histc');
xlabel('\Delta Distance from Trail Between Sniffs (mm)');
ylabel('Probability of Turn');