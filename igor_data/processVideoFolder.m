function vids = processVideoFolder(folder_name, tracker, varargin)
% function vids = processVideoFolder(folder_name, tracker)
%
% Processes a folder of video files using the type of tracker object that is given as TRACKER in the
% input.  VARARGIN is used to specify whether or not to save the files - you can call this to load a
% group of objects in a folder easily, and not saving speeds that operation up.  Obviously, the default
% is to have the the processed tracker objects save themselves.

if ~isempty(varargin)
    saveFlag = varargin{1};
else 
    saveFlag = 1;
end

global VIDEO_ROOT;
folder_path = [VIDEO_ROOT filesep folder_name filesep];
ind_name = [folder_path 'tracking_times.txt'];
if ~exist(ind_name, 'file')
    disp('There is no index file for the videos');
else
    starts = []; ends = [];
    fid = fopen(ind_name, 'r');
    res = textscan(fid, '%s [%d %d]');
    fnames = res{1};
    starts = res{2};
    ends = res{3};
    
    s = matlabpool('size');
    if s==0 %check if parallel toolbox is running.  If not, just do regular for loop
        for ii = 1:length(fnames)
            vids(ii) = trackAndSave(tracker, folder_path, fnames{ii}, starts, ends, ii, saveFlag);
        end
    else
        parfor ii = 1:length(fnames)
            vids(ii) = trackAndSave(tracker, folder_path, fnames{ii}, starts, ends, ii, saveFlag);
        end
    end
end

% ------------------------ Helper function ---------------------
function vid = trackAndSave(tracker, folder_path, fname, starts, ends, ii, saveFlag)
    if (starts(ii)==0 && ends(ii)==0) %the tag to just read the whole video. Leaving blank destroys file reading
        vid = tracker([folder_path fname], []);
    else
        vid = tracker([folder_path fname], [], [starts(ii), ends(ii)]); %consider the given range
    end
    vid.mousePosition([]); %compute position in whole movie, as defined on opening
    if saveFlag
        vid.save(); %save tracking to mat file
    end