function [i1,i2, p1, p2] = intersectPosition(v1, v2)
% function [i1,i2] = intersectPosition(v1, v2)
% 
% Returns the indices of intersection of the two vectors, or the closest point to intersection
% The inputs, v1 and v2, must be n x 2 matrices.

% Compute the interpoint distance matrix, then just return the minimum indices and values
dm = ipdm(v1,v2);

[d1,i1] = min(dm,[],2); 
[d2,i2] = min(dm,[],1);

p1 = v1(i1,:);
p2 = v2(i2,:);



