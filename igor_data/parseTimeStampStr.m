function ts = parseTimeStampStr(timeStr)
% function parseTimeStampStr(tsStr)
% 
% Parses a string of the the format 'TIME=xxxxxx.xxxx' and returns the
% floating point value of the xxxx.xxxx

remain = timeStr; toks = {}; ii=1;
while ~isempty(remain)
    [toks{ii} remain] = strtok(remain, '=');
    ii = ii+1;
end
ts = str2double(toks{end});