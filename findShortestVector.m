function closestP = findShortestVector(points, verts)
% function [dists, vects] = findShortestVector(points, verts)
%
% The point of this is to find the shortest vector to a shape defined by vertices on a grid.
% This works in a specific case - where it is given a short set of points (verts) constituing the closest points
% each input point.  If there are father verts are given, then it will give nonsense out. 

closestP = verts(:,:,1);
points3 = repmat(points, [1 1 size(verts,3)]); %copies by how many points there are 
less = verts-points3; 
close = abs(less) <= 1;

for ii=1:2 %x and y dimensions
    nclose = sum(close(:,ii,:),3);
    between = nclose >= 2;
    closestP(between, ii) = points(between,ii);
end
    
