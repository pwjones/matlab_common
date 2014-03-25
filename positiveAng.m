function a = positiveAng(a)
% function a = positiveAng(a)
%
% utility function that makes negative angles positive in the range of 0-2pi
np = -floor(a./(2*pi));
pos = np>0; %selection logical array based on those indices that were positive (originally negative values)
a(pos) = a(pos)+(np(pos).*(2*pi));
a = mod(a, 2*pi);
