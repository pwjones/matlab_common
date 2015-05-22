function [paraRel, absAng] = convertMouseHeadingAngles(orthoRel, ortho, pos)
% function para_angles = convertMouseHeadingAngles(ortho_angles, pos)
%
%
% The purpose of this is to convert from Orthogonal angles to parallel 
% angles.
absAng = orthoRel + ortho;
isright = pos > 0;
para = ortho;
para(isright) = ortho(isright) + pi/2; 
para(~isright) = ortho(~isright) + 3*pi/2;
paraRel = circ_dist(absAng, para);