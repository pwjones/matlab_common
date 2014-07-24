% Code to test the similarity (overlap/position differe) of blobs in adjacent frames

%vid = vids(8);

nf = vid.nFrames;
dists = NaN*zeros(nf*(size(vid.areas,2)^2),1); overlap = NaN*zeros(nf*(size(vid.areas,2)^2),1);
size_diff = NaN*zeros(nf*(size(vid.areas,2)^2),1);
disti = 1;
for ii = 1:(nf-1)
   area = vid.areas(ii,:);
   area1 = vid.areas(ii+1,:);
   for jj = 1:vid.nblobs(ii)
       % getting the distance between centroids
        nb = vid.nblobs(ii+1);
        if (disti+nb) > nf
            break;
        end
        cents = [area1.Centroid];
        cents = cents(1:(vid.nblobs(ii+1)*2));
        cents = reshape(cents, 2,[])'; 
        c = repmat(area(jj).Centroid, vid.nblobs(ii+1),1);
        dist = sqrt(sum((cents-c).^2,2));
        dists(disti:(disti+nb-1)) = dist(:);
        % getting the percentage overlap between blobs
        ol = regionOverlap(area(jj),area1, {'PixelIdxList'});
        overlap(disti:(disti+nb-1)) = ol(1:nb);
        sd = abs([area1(1:nb).Area] - area(jj).Area)/area(jj).Area;
        size_diff(disti:(disti+nb-1)) = sd(:);
        disti = disti+nb;
   end
end

figure;
hist(dists, 500);
xlabel('Distance between Centroids (px)');
figure;
hist(overlap, 500);
xlabel('Overlap of regions (proportion)');
figure;
hist(size_diff, 500);
xlabel('Size difference between blobs (px)');

%% Let's do this again, but now just look at the ones that are identified by the algorithm as "nose".
 % This may be a bit circular, but the distributions within the ranges should be useful.
 
nf = vid.nFrames;
dists = NaN*zeros(nf*size(vid.areas,2),1); overlap = NaN*zeros(nf*size(vid.areas,2),1);
size_diff = NaN*zeros(nf*size(vid.areas,2),1);
disti = 1;
for ii = 1:(nf-1)
   area = vid.areas(ii,:);
   area1 = vid.areas(ii+1,:);
   nb = vid.noseblob(ii); nb1 = vid.noseblob(ii+1);
   if ~isnan(nb) && ~isnan(nb1)
       dist = area1(vid.noseblob(ii+1)).Centroid - area(vid.noseblob(ii)).Centroid;
       dist = sqrt(sum(dist.^2));
       dists(ii) = dist;
       overlap(ii) = regionOverlap(area(nb),area1(nb1), {'PixelIdxList'});
       size_diff(ii) = abs(area1(nb1).Area - area(nb).Area)/area(nb).Area;
   end
end

figure;
hist(dists, 100);
xlabel('Distance between Centroids (px)');
figure;
hist(overlap, 100);
xlabel('Overlap of regions (proportion)');
figure;
hist(size_diff, 100);
xlabel('Size difference between blobs (px)');


%%

tracker = vids(11);

paths = cat(1, tracker.paths(1).PixelList, tracker.paths(2).PixelList);
bodyCOM = tracker.bodyCOM;
dm = ipdm(paths,bodyCOM);
min_dist = min(dm);
close = min_dist < 30;
missing = isnan(tracker.nosePos(:,1));
inds = close & missing;



