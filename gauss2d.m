function weights = gauss2d(center, sigma2, pos_x, pos_y)
% function weights = gauss2d(center, sigma2, pos_x, pos_y)
% Returns the values of a 2D gaussian, where the center (X,Y) coordinates are given in 
% CENTER, and the variance of both dimensions are given by SIGMA2.  The positions for
% which to return values are specified in POS_X and POS_Y.

center_x = center(1);
center_y = center(2);
weights = 1/(2*pi* sigma2).*exp(-1*((pos_x-center_x).^2+(pos_y-center_y).^2)./(2* sigma2));