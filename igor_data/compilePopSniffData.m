function sniffData = compilePopSniffData()

% Load all of the processed data from the Spring-Fall 2014 sniffing dataset
global SNIFF_ROOT;
base_folder = SNIFF_ROOT;
save_flag = 1;
mouse_names = {'19439', '21413', '971', '1080', '1527'};
%mouse_names = {'19439'};
shiftSniffs = 0;
folders = {''};
nMice = length(mouse_names);
sniffData = [];
for ii = 1:nMice
    saved_file = [base_folder filesep mouse_names{ii} '.mat'];
    % Each file is for a mouse
    if exist(saved_file, 'file')
        disp(['Loading: ' saved_file]);
        load(saved_file);
        % The variable 'mouseData' is created in the variable space
        sniffData = extractMouseSniffData(mouseData, sniffData, shiftSniffs);
    end
end



% ------------------------------------------------
function data = extractMouseSniffData(exp, data, shiftSniffs)
% function data = extractMouseSniffData(exp, data)
%
% Called once per each mouse (exp structure).  Then appends to whatever
% data is present in the 'data' input structure.  

% Parameters
thresh_dist = 20; %mm
thresh_vel = 40; %mm/sec
mm_conv = .862; %mm/px linear
traj_wind = -10:15; % #of samples around the turning to retrieve - 50fps, 20ms per frame
turn_wind = .08; %s, time after a sniff to look for turns
nboot = 1000;

% data structures for the turning probabilities
dist_diffs = []; bturn = []; turnDirections = []; dirToTrail = []; mouseHeadingFromOrtho = []; 
allPreTurnHeadings = []; allTurnTrajectories = []; allTurnDirs = []; 
allPostTurnHeadings = []; allPreTurnDistDiff = []; turnLag = [];
allSniffDirToTrail = []; allSniffHeadingsFromOrtho = []; allSniffISI = []; allTurnSniffPos = [];
newSeg = []; allTurnPos = []; allSniffPos = []; allPostTurnDirToTrail = []; allPreTurnDirToTrail = [];
allPostTurnPos = []; allPrePos = [];

for ii = 1:length(exp.vids)
    exp.vids(ii).makePathsSkel();
    % Velocities
    noseVel = exp.vids(ii).noseVel * mm_conv * exp.vids(ii).frameRate; %get the nose/body velocities
    bodyVel = exp.vids(ii).bodyVel(:,1) * mm_conv * exp.vids(ii).frameRate;
    
    dt = exp.vids(ii).times(2) - exp.vids(ii).times(1);
    
    % Let's gather the sniff information
    sniffFrames = exp.resp(ii).sniffFrames(exp.resp(ii).vidSniffs);
    sniffTimes = exp.resp(ii).sniffTimes(exp.resp(ii).vidSniffs);
    if shiftSniffs
        shifts = randi([0 ceil(mean(diff(sniffFrames)))], [nboot, 1]);
    end
    
    %The nose trajectories and following segments
    odist= exp.vids(ii).orthogonalDistFromTrail(1:exp.vids(ii).nFrames, 1);
    [headingFromOrtho, ~, dtt] = exp.vids(ii).headingFromMotion_TrailRelative(1:exp.vids(ii).nFrames,1);
    odist = double(mm_conv * odist);
    fsegs = exp.vids(ii).getFollowingSegments([], 1, thresh_dist);
    
    % Get the turning information - these are turns that happen within the
    % thresh_dist from the trail, so it's already selected for following
    [turningInds, dirs, nose_traj] = exp.vids(ii).findFollowingTurns([], 1, thresh_dist, traj_wind);
    % Need to just exclude turns that are too early or late for proximal
    % events to be out of range
    sel = (turningInds + traj_wind(1)) > 0  & (turningInds - traj_wind(1) <= exp.vids(ii).nFrames);
    turningInds = turningInds(sel); dirs = dirs(sel); nose_traj = nose_traj(:,sel);
    
    % Finding the distance change between previous 2 sniffs
    dd = NaN*zeros(length(turningInds), 3);
    sniffPos = NaN*ones(length(turningInds), 4); %4 previous sniffs
    lag = NaN*zeros(length(turningInds),1);
    sniffi = NaN*ones(length(turningInds),1);
    
    for jj=1:length(turningInds)
        preTS = find(sniffFrames < turningInds(jj), 4, 'last');
        if length(preTS) == 4
            dd(jj,:) = diff(abs(odist(sniffFrames(preTS))));
            lag(jj) = turningInds(jj) - sniffFrames(preTS(end));
            sniffi(jj) = sniffFrames(preTS(end));
            sniffPos(jj,:) = odist(sniffFrames(preTS));
            [sniffHeadingFromOrtho(jj,:), ~, sniffDirToTrail(jj,:)] = exp.vids(ii).headingFromMotion_TrailRelative(sniffFrames(preTS),1);
        end
    end
    sel = ~isnan(sniffi); % select and trim possible NaNs.
    turningInds = turningInds(sel); dirs = dirs(sel); nose_traj = nose_traj(:,sel);
    dd = dd(sel,:); lag = lag(sel); sniffi = sniffi(sel); sniffPos = sniffPos(sel,:);
    sniffDirToTrail = sniffDirToTrail(sel, :);
    sniffHeadingFromOrtho = sniffHeadingFromOrtho(sel, :);
    
    tpos = exp.vids(ii).orthogonalDistFromTrail(turningInds,1); % position at turning frame
    [thead, ~, toTrail] = exp.vids(ii).headingFromMotion_TrailRelative(turningInds-4, 1); %heading before turn
    prePos = exp.vids(ii).orthogonalDistFromTrail(turningInds-4, 1);
    [theadPost, ~, toTrailPost] = exp.vids(ii).headingFromMotion_TrailRelative(turningInds+4, 1); %heading after turn
    postPos = exp.vids(ii).orthogonalDistFromTrail(turningInds+4, 1);
    %sniffhead = exp.vids(ii).headingFromMotion_TrailRelative(sniffi, 1); %heading at the sniff preceding
    
    tpos = double(mm_conv.*tpos);
    allTurnTrajectories = [allTurnTrajectories, nose_traj];
    allPreTurnHeadings = [allPreTurnHeadings; thead];
    allPreTurnDirToTrail = [allPreTurnDirToTrail; toTrail];
    allPrePos = [allPrePos; prePos];
    allSniffDirToTrail = [allSniffDirToTrail; sniffDirToTrail];
    allSniffHeadingsFromOrtho = [allSniffHeadingsFromOrtho; sniffHeadingFromOrtho];
    allPostTurnHeadings = [allPostTurnHeadings; theadPost];
    allPostTurnDirToTrail = [allPostTurnDirToTrail; toTrailPost];
    allPostTurnPos = [allPostTurnPos; postPos];
    allTurnDirs = [allTurnDirs; dirs(:)];
    allPreTurnDistDiff = [allPreTurnDistDiff; dd];
    allTurnSniffPos = [allTurnSniffPos; sniffPos];
    allTurnPos = [allTurnPos; tpos];
    
    
    % Its best if we do this following segment by segment, because if
    % we anneal them together and try to look at concentration differences
    % between sniffs then we'll get unrelated sniff concentration
    % differences.
    followi = [];
    for jj = 1:size(fsegs,1)
        range = fsegs(jj,1):fsegs(jj,2); %range of frames
        [turns, ~,ti] = intersect(range, turningInds); % selecting within range items
        [sniffs,~,sniffind] = intersect(range, sniffFrames);
        turnDirs = dirs(ti); %turning directions
        sniffT = sniffTimes(sniffind);
        sniffISI = diff(sniffT); %this is the sniffing during following
        
        for kk=1:(length(sniffs))
            sniffFrameRange = sniffind(kk)-3:sniffind(kk); %indices into the sniffFrames vector
            sniffFrameRange(sniffFrameRange < 1) = 1; %make sure it doesn't fall below 1
            pos = odist(sniffFrames(sniffFrameRange)); 
            allSniffPos = [allSniffPos; pos'];
            diff_dist = diff(abs(pos)); diff_dist(diff_dist==0) = NaN; %distance between sniff pos then NaN any 0 from missing sniff positions
            dist_diffs = [dist_diffs; diff_dist'];
            
            %dist_diffs = [dist_diffs; abs(odist(sniffs(kk+1))) - abs(odist(sniffs(kk)))];
            dirToTrail = [dirToTrail; dtt(sniffFrameRange)']; %the mouse angle relative to the vector to trail
            mouseHeadingFromOrtho = [mouseHeadingFromOrtho; headingFromOrtho(sniffFrameRange)'];
            frame_wind = sniffs(kk):(sniffs(kk) + floor(turn_wind/dt));
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
            % away from the trail. NOTE THAT THIS IS DIFFERENT THAN THE PROCESSING
            % DONE FOR THE NON-SNIFF POPULATION DATA - THAT IS KEPT AS LEFT/RIGHT.
            % ALSO THIS IS USING THE sniff position rather than the actual turn position
            %turnDirections = [turnDirections; td*odist(sniffs(kk))];
            turnDirections = [turnDirections; td];
            if (kk == 1) newSeg = [newSeg; 1]; else newSeg = [newSeg; 0]; end
        end
        %dist_diffs = [dist_diffs; NaN, NaN, NaN]; %want to pad one row of NaNs between segments
        allSniffISI = [allSniffISI; sniffISI(:)];
        followi = [followi, range];
    end
    
end

if isempty(data)
    % Assembled by looking at each following section individually.  
    data.dist_diffs = []; 
    data.segmentBegin = [];
    data.bturn = []; 
    data.turnDirections = []; 
    data.dirToTrail = [];
    data.mouseHeadingFromOrtho = [];
    data.followingSniffISI = [];
    data.sniffPos = [];
    % Assembled by getting the turn info first, basing events around the
    % turns.
    data.turnTrig_preTurnHeadings = [];
    data.turnTrig_preTurnDirToTrail = [];
    data.turnTrig_preTurnPos = [];
    data.turnTrig_turnTrajectories = []; 
    data.turnTrig_turnDirs = []; 
    data.turnTrig_postTurnHeadings = []; 
    data.turnTrig_postTurnDirToTrail = [];
    data.turnTrig_postTurnPos = [];
    data.turnTrig_preTurnDistDiff = []; 
    data.turnTrig_turnLag = [];
    data.turnTrig_sniffHeadings = [];
    data.turnTrig_sniffDirToTrail = [];
    data.turnTrig_sniffPos = [];
    data.turnTrig_turnPos = [];
end
% Concatenate results - this seems like the best solution, as analysis can
% treat the results of multiple extractions as single data arrays within the structure.
data.dist_diffs = cat(1, data.dist_diffs, dist_diffs);
data.segmentBegin = cat(1, data.segmentBegin, newSeg);
data.bturn = cat(1, data.bturn, bturn);
data.turnDirections = cat(1, data.turnDirections, turnDirections);
data.dirToTrail = cat(1,data.dirToTrail, dirToTrail);
data.mouseHeadingFromOrtho = cat(1,data.mouseHeadingFromOrtho, mouseHeadingFromOrtho);
data.sniffPos = cat(1,data.sniffPos, allSniffPos);
data.turnTrig_preTurnHeadings = cat(1,data.turnTrig_preTurnHeadings, allPreTurnHeadings);
data.turnTrig_preTurnDirToTrail = cat(1,data.turnTrig_preTurnDirToTrail, allPreTurnDirToTrail);
data.turnTrig_preTurnPos = cat(1, data.turnTrig_preTurnPos, allPrePos);
data.turnTrig_turnTrajectories = cat(2,data.turnTrig_turnTrajectories, allTurnTrajectories);
data.turnTrig_turnDirs = cat(1,data.turnTrig_turnDirs, allTurnDirs);
data.turnTrig_postTurnHeadings = cat(1,data.turnTrig_postTurnHeadings, allPostTurnHeadings);
data.turnTrig_postTurnDirToTrail = cat(1, data.turnTrig_postTurnDirToTrail, allPostTurnDirToTrail);
data.turnTrig_postTurnPos = cat(1,data.turnTrig_postTurnPos, allPostTurnPos);
data.turnTrig_preTurnDistDiff = cat(1,data.turnTrig_preTurnDistDiff, allPreTurnDistDiff);
data.turnTrig_turnLag = cat(1,data.turnTrig_turnLag, turnLag);
data.turnTrig_sniffHeadings = cat(1,data.turnTrig_sniffHeadings, allSniffHeadingsFromOrtho);
data.turnTrig_sniffDirToTrail = cat(1,data.turnTrig_sniffDirToTrail, allSniffDirToTrail);
data.followingSniffISI = cat(1, data.followingSniffISI, allSniffISI);
data.turnTrig_sniffPos = cat(1, data.turnTrig_sniffPos, allTurnSniffPos);
data.turnTrig_turnPos = cat(1, data.turnTrig_turnPos, allTurnPos);
data.traj_wind = traj_wind;