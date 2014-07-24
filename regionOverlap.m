function overlap = regionOverlap(area1, areas2, regprops)
% function regionSimilarity(area1, areas2, regprops)
%
% Trying to come up with a similarity score for the two areas

% regprops = {'Area', 'Centroid', 'BoundingBox','MajorAxisLength','MinorAxisLength','Orientation','Extrema','PixelIdxList','PixelList','Perimeter'};
dbg = 0;
if dbg
    im = zeros(1024,1280);
end
overlap = zeros(length(areas2),1);
if sum(strcmp('PixelIdxList', regprops))
    px1 = area1.PixelIdxList;
    if dbg im(px1) = im(px1)+1; end
    for ii = 1:length(areas2)
        px2 = areas2(ii).PixelIdxList;
        if dbg im(px2) = im(px2)+1; end
        overlap(ii) = length(intersect(px1, px2)) / ((length(px1)+length(px2))/2);
    end
end

if dbg
    figure;
    imshow(im);
end
