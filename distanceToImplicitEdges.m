function [dist, vect] = distanceToImplicitEdges(pt, verts)
% pt should be an N dim pt (that is a floating point).
% Verts should be an M x N array of M whole number points, each in N
% dimensions.

new_verts = verts;
for i = 1:length(pt) %deal with each dim separately
    c = ceil(pt(i));
    f = floor(pt(i));
    for j = 1:size(verts,1)
        