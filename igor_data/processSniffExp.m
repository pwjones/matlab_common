function exp = processSniffExp(expFilename, vidBasename)
% function exp = processSniffExp(expFilename, vidBasename)
%
% Process a set of data files that includes a single EXPFILENAME.h5 file
% that has been exported using ExportHDF5() from Igor that contains a set
% of sweeps and the VIDBASENAME*.avi movie files that correspond to those 
% sweeps.  Assumes that the video captured was triggered by a trigger signal present
% in the h5 file. Detects events in the physiology and tracks the animal movement using 
% the videos. MUST CURRENTLY BE RUN FROM THE DIRECTORY WHERE THE FILES ARE
% LOCATED. 
% OUTPUT: a single structure for the processed experimental data, EXP.

dbg=0; %debug flag
expFN = [expFilename, '.h5'];
% Setting initial vars
respSignal =  'Thermocouple';
%tracker = @MouseTrackerUnder2;
tracker = @MouseTrackerKF
% Filter lengths, now in time units
median_filt_len = 3e-3; %median filter time in sec
hp_len = 40e-3; % high-pass filter the signal to remove baseline fluctuations, sec 
%median_filt_len = 30; % number of samples to median filter over - denoise the signal a little
%hp_len = 400; % high-pass filter the signal to remove baseline fluctuations, # of samples 

%% %%%%%%%%%%%%%%%%%%%%%%%%% Load and process the videos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%vidBasename = '11963';

files = dir([vidBasename '*.avi']); %finds the files that start with that name
% Now we need to see if there is a file specifying the ranges of the movies
% to analyze
starts = []; ends = [];
time_file = dir('tracking_times.txt');
if ~isempty(time_file)
    fid = fopen(time_file.name, 'r');
    res = textscan(fid, '%s [%d %d]');
    fnames = res{1};
    starts = res{2};
    ends = res{3};
    % select the files that start with vidBasename
    inc = zeros(length(fnames),1);
    for ii=1:length(fnames)
        tok = strtok(fnames{ii}, '_');
        if strcmp(tok, vidBasename)
            inc(ii) = 1;
        end
    end
    inc = logical(inc);
    fnames = fnames(inc); starts = starts(inc); ends = ends(inc);
end
for ii = 1:length(fnames)
    fname = fnames{ii};
    if (starts(ii)==0 && ends(ii)==0) %the tag to just read the whole video. Leaving blank destroys file reading
        vid = tracker(fname, []);
    else
        vid = tracker(fname, [], [starts(ii), ends(ii)]); %consider the given range
    end
    vid.mousePosition([]); %compute position in whole movie, as defined on opening
    vid.save(); %save tracking to mat file
    vids(ii) = vid;
end
disp('Checking for missing frames in videos'); %doing exactly that
for ii = 1:length(vids)
    missing = vids(ii).isMissingFrame();
    if missing
        disp(sprintf('%s is missing a frame.', vids(ii).videoFN));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Camera trigger pulse segment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
camSignal = 'Camera';
%vids= [];
camTrig = readIgorH5file(expFN, camSignal, 0, 0);
% loop to find the frame indices - more explicitly the first sample of each high pulse in the trigger signal 
camTrig = processCameraPulses(camTrig, vids, dbg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LED pulses signaling the frame timing %%%%%%%%%%%%%%%%%%%%%%%%%%%%
ledStr = 'LED';
led = readIgorH5file(expFN, ledStr, 0, 0);
% We need to do something else to verify timing.


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%Load the respiration data from a thermocouple %%%%%%%%%%%%%%%%%%%%
resp = readIgorH5file(expFN, respSignal, median_filt_len, hp_len);
resp = computeSniffTimes(resp);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Select particular trials if necessary %%%%%%%%%%%%%%%%%%%%%%%%%%
trials_file = dir('trials.txt');
if ~isempty(trials_file)
    fid = fopen('trials.txt');
    finds = fscanf(fid, '%d');
    vids = vids(finds);
    resp = resp(finds);
    camTrig = camTrig(finds);
end
%% Define returned structure
exp.vids = vids;
exp.resp = resp;
exp.camTrig = camTrig;

exp = findSniffFrames(exp);
