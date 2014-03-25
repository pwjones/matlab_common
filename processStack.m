function [bright_vect, area] = processStack(base_filename)

pb=1;

f = dir([base_filename '*']);
bright_vect = zeros(length(f), 1);
for i = 1:length(f)
    im = imread(f(i).name);
    if(i==1)
        figure; ih = imshow(im);
        area = impoly(gca);
        im_mask = area.createMask();
    end
    gim = im(:,:,1);
    px = gim(im_mask);
    bright_vect(i) = mean(px(:));
end

if pb
    figure; hold on;
    plot(bright_vect);
    xlabel('Sample');
    ylabel('Average Brightness');
    
    figure;
    im = imread(f(1).name);
    im2_layer = im(:,:,1);
    im2_layer(im_mask)=255;
    im2 = cat(3,im2_layer, im(:,:,2), im(:,:,3));
    imshow(im2);
    
    figure;
    im = imread(f(end).name);
    im2_layer = im(:,:,1);
    im2_layer(im_mask)=255;
    im2 = cat(3,im2_layer, im(:,:,2), im(:,:,3));
    imshow(im2);
end


        
        