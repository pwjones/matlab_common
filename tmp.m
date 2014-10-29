% temporary scratch code
overlay2 = this.trackingStruct.binMovie(:,:,(framei-1)*this.trackingStruct.movieSubsample +3);
g = f;
r = cast(overlay2,'uint8')*255;
g(:,:,1) = g(:,:,1) + r;
imshow(g);