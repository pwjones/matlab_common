function contained = isContained(points, verts)
% function isContained(testPoints, verts)
% 
% Trying to see if testPoints are contained within the verts, which are the 4 closest points to each 
% test point.  
% n x 2 for testPoints
% n x 2 x 4 for verts

points3 = repmat(points, [1 1 size(verts,3)]); %copies by how many points there are 
less = verts-points3; 
close = abs(less(:,1,:)) <= 1 & abs(less(:,2,:)) <=1;
nclose = squeeze(sum(close,3));
contained = nclose >= 4;

% less = squeeze(reshape(less, size(less, 1), [], 1));
% maxdist = nanmax(abs(less),[], 2);
% contained = maxdist <= 1;

