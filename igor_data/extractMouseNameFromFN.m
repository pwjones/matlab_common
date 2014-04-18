function mouseName = extractMouseNameFromFN(filename)

remain = filename;
while ~isempty(remain)
    [str remain] = strtok(remain,filesep);
end
%str should be the file without folders
if ~isempty(str)
    mouseName = strtok(str, '_');
else
    disp('extractMouseNameFromFN: Couldnt find the mouse name in given string');
    mouseName = '';
end