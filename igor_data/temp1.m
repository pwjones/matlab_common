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
end

for ii = 1:length(fnames)
    fname = fnames{ii};
    if (starts(ii)==0 && ends(ii)==0) %the tag to just read the whole video - leaving blank destroys file reading
        vid = MouseTrackerUnder2(fname, []);
    else
        vid = MouseTrackerUnder2(fname, [], [starts(ii), ends(ii)]); %consider the given range
    end
    vid.mousePosition([]); %compute position in whole movie, as defined on opening
    vid.save(); %save tracking to mat file
    vids(ii) = vid;
end

% Let's do some processing 

% Let's get the frames when the animal sniffed
for ii=1:length(exp.resp)
    frame_t = exp.camTrig(ii).time(exp.camTrig(ii).frameInds);
    sniff_times = exp.resp(ii).time(exp.resp(ii).sniffVect);
    sniff_frames = zeros(length(sniff_times),1);
    for jj=1:length(sniff_times)
        sniff_frames(jj) = find(frame_t >= sniff_times(jj), 1, 'first');
    end
    exp.resp(ii).sniff_times = sniff_times;
    exp.resp(ii).sniff_frames = sniff_frames;
    exp.resp(ii).sniff_pos = exp.vid(ii).nosePos(sniff_frames,:);
end

%%
gt = .09;
bm = false(size(mf));
for ii = 1:size(mf,3) 
    bw = im2bw(mf(:,:,ii),gt);
    bw = imerode(bw, strel('square',3)); 
    bm(:,:,ii) = bw;
end
%%

cc = bwconncomp(bm);
stats = regionprops(cc, 'Area');
idx = find([stats.Area] > 20);
bm = ismember(labelmatrix(cc), idx);
