function cm = blackGradColormap(color, n)
% function cm = blackGradColormap(color, n)
%
% Generates a colormap where the first entry is black and the others are evenly spaced to a given color.
% Useful for visualizing different levels of a thing.

inc = color/(n-1);
cm = zeros(n,3);
for ii = 1:(n-1)
    cm(ii+1,:) = inc*ii;
end
% error checking
toohigh = cm > 1;
cm(toohigh) = 1;
low = cm < 0;
cm(low) = 0;

