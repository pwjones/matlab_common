function rot_dir = rotationDirection(ref_vect, test_vect)
% function rot_dir = rotationDirection(ref_vect, test_vect)
% 
% Given a reference vector (ref_vect), which can be a vector of 
% directions, from 0-2pi radians, give the handedness of the rotation to
% the test vector (test_vect).  Return ROT_DIR which is 1 for right handed
% rotation, -1 for left handed rotation, and 0 for no rotation.

% make all of the rotations positive, 0-2pi
ref_vect = mod(ref_vect,2*pi);
rot_norm = pi - ref_vect;
ref_rot = ref_vect + rot_norm; %is this always pi? should be.

%rotate the test the same way
test_vect = mod(test_vect, 2*pi);
test_rot = mod(test_vect+rot_norm,2*pi);

%assign the output
rot = test_rot < pi;
rot_dir = double(rot);
rot_dir(~rot) = -1;
same = (test_rot == pi) | (test_rot == 0);
rot_dir(same) = 0;

