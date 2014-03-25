function overlap = regionOverlap(area1, areas2, regprops)
% function regionSimilarity(area1, areas2, regprops)
%
% Trying to come up with a similarity score for the two areas

% regprops = {'Area', 'Centroid', 'BoundingBox','MajorAxisLength','MinorAxisLength','Orientation','Extrema','PixelIdxList','PixelList','Perimeter'};


overlap = zeros(length(areas2),1);
if sum(strcmp('PixelIdxList', regprops))
    px1 = area1.PixelIdxList;
    for ii = 1:length(areas2)
        px2 = areas2(ii).PixelIdxList;
        overlap(ii) = length(intersect(px1, px2)) / ((length(px1)+length(px2))/2);
    end
end