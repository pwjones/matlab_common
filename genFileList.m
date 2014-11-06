function out_str = genFileList(folder, varargin)
% Just a function to list the names of files in a directory.
% Useful for copying and pasting

global VIDEO_ROOT;
sort_by_timestamp = 1;

if ispc
    base_folder = 'I:\\odor_tracking'
else
    base_folder = VIDEO_ROOT;
end
if nargin>1
    match_str = varargin{1};
else
    match_str = [];
end
ts = {};
if sort_by_timestamp
    match_str = ['*' match_str '*'];
    files = dir([base_folder, filesep, folder, filesep, match_str]);
    out_str = [];
    delim = '_';
    for i=1:length(files)
        name = files(i).name;
        [~, name] = strtok(name, delim);
        ts{i} = strtok(name, '._');
    end
    [~, sorti] = sort(ts);
    files = files(sorti);
end
for i=1:length(files)
    temp_str = [files(i).name ' []'];
    disp(temp_str);
    %out_str = [out_str; temp_str];
end
