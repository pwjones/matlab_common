function newA = increaseMovContrast(imarray)
%function newA = increaseMovContrast(imarray)
% 
% Stretches the image histogram to be the limits of the entire movie.
% Assumes 8 bit images 
mov_max = single(max(imarray(:)));
mov_min = single(min(imarray(:)));

newA = zeros(size(imarray), 'uint8');
for ii = 1:size(imarray,3)
    %mov_max = single(max(max(imarray(:,:,ii))));
    newA(:,:,ii) = imadjust(imarray(:,:,ii), [mov_min/255 mov_max/255], [0 1]);
end
%imarray = uint8(newA*255);
