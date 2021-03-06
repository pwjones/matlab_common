function vid = trackAndSave(tracker, folder_path, fname, starts, ends, ii, saveFlag)
% function vid = trackAndSave(tracker, folder_path, fname, starts, ends, ii, saveFlag)
%
% ------------------------ Helper function to processVideoFolder ---------------------
disp(['In trackAndSave: ' folder_path fname]);
if isempty(starts)
    starts = zeros(ii,1);
    ends = zeros(ii,1);
end
if (starts(ii)==0 && ends(ii)==0) %the tag to just read the whole video. Leaving blank destroys file reading
    disp('Calling without limits');
    vid = tracker([folder_path fname], []);
else
    disp('Calling with limits');
    vid = tracker([folder_path fname], [], [starts(ii), ends(ii)]); %consider the given range
end
if (vid.mm_conv == 1.16) % I had the wrong unit conversion for awhile, this will help fix that for already saved files
    vid.mm_conv = .862;
    saveFlag = 1;
end
vid.mousePosition([]); %compute position in whole movie, as defined on opening
if saveFlag
    vid.save(); %save tracking to mat file
end