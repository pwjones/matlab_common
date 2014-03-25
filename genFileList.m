function out_str = genFileList(folder)
% Just a function to list the names of files in a directory.
% Useful for copying and pasting

global VIDEO_ROOT;

if ispc
    base_folder = 'I:\\odor_tracking'
else
    base_folder = VIDEO_ROOT;
end

files = dir([base_folder, filesep, folder]);
out_str = [];
for i=1:length(files)
    temp_str = [files(i).name ' []'];
    disp(temp_str);
    %out_str = [out_str; temp_str];
end
