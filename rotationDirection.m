function [rot_dir, test_rot]  = rotationDirection(ref_vect, test_vect)
% function rot_dir = rotationDirection(ref_vect, test_vect)
% 
% Given a reference vector (ref_vect), which can be a vector of 
% directions, from 0-2pi radians, give the handedness of the rotation to
% the test vector (test_vect), which is also given by values from 0-2pi.  
% These are vectors in the array sense, not in the magnitude/direction sense.
% Both inputs are directions.  
% Returns ROT_DIR which is 1 for right handed
% rotation, -1 for left handed rotation, and 0 for no rotation.
%
% Test_rot is the angular difference between the two angles, with 0 being
% the equality and ranging from -pi:pi


ref_vect = mod(ref_vect,2*pi); % make all of the rotations positive, 0-2pi
ref_norm = pi - ref_vect; %now shift so that the center is back at 0; -pi:pi range
%pi_check = ref_vect + ref_norm; %is this always pi? should be.

%rotate the test to be 0:2pi
test_vect = mod(test_vect, 2*pi);
% sum them together: >pi is counter-clockwise, <pi is clockwise
test_rot = mod(test_vect+ref_norm,2*pi);

%assign the output
rot = test_rot < pi; %rightward rotation (clockwise)
rot_dir = double(rot);
rot_dir(~rot) = -1; %mark leftward (counter-clockwise) with negative
same = (test_rot == pi) | (test_rot == 0);
rot_dir(same) = 0;
test_rot = test_rot-pi;

