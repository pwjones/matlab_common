function [allMinis, allMiniSlopes, allMiniHeights] = returnAllMinis(handles)
% function [allMinis, allMiniSlopes] = returnAllMinis(handles)
% 
% Purpose of the function is to give structures of all minis, current and
% previous
if(~isempty(handles.prev))
    prevNum = sum([handles.prev.numMinis]);
else
    prevNum = 0;
end
numMinis = handles.numMinis + prevNum; % total mini count
allMinis = NaN * zeros(size(handles.minis, 1), numMinis);
allMiniSlopes = NaN * zeros(size(allMinis,1)-1, numMinis);
allMiniHeights = NaN * zeros(numMinis, 1);
mini_n = 0;
for i=1:length(handles.prev)
    nums = mini_n + (1:handles.prev(i).numMinis);
    allMinis(:,nums) = handles.prev(i).minis(:, 1:handles.prev(i).numMinis);
    allMiniSlopes(:,nums) = handles.prev(i).miniSlopes(:, 1:handles.prev(i).numMinis);
    allMiniHeights(nums) = handles.prev(i).miniHeights(1:handles.prev(i).numMinis);
    mini_n = nums(end);
end
nums = mini_n + (1:handles.numMinis);
allMinis(:,nums) = handles.minis(:,1:handles.numMinis);
allMiniSlopes(:,nums) = handles.miniSlopes(:, 1:handles.numMinis);
allMiniHeights(nums) = handles.miniHeights(1:handles.numMinis);