function vert = pixelVerticesFromPoints(pixelList)

px = pixelList; %n x 2
vert = NaN*zeros(size(px,1)*4, 2); %upper bound is 4 times the px
vi = 1;
for jj = 1:size(px,1)
    pxt = px(jj,:);
    vert(vi,:)   = [pxt(1)-.5, pxt(2)-.5]; %left, top
    vert(vi+1,:) = [pxt(1)+.5, pxt(2)-.5]; %right, top
    vert(vi+2,:) = [pxt(1)-.5, pxt(2)+.5]; %left, bottom
    vert(vi+3,:) = [pxt(1)+.5, pxt(2)+.5]; %right, bottom
    vi = vi+4;
end
vert = unique(vert, 'rows', 'R2012a');
