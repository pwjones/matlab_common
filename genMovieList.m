%generate a text output 

function out_str = genMovieList(folder)

if ispc
    base_folder = 'I:\\odor_tracking'
else
    base_folder = '/Users/pwjones/Movies/mouse_training';
end

files = dir([base_folder, filesep, folder]);
out_str = [];
for i=1:length(files)
    temp_str = [folder filesep files(i).name ' []'];
    disp(temp_str);
    %out_str = [out_str; temp_str];
end
