function X = bootrsmp(n, y, q, L)
% function bootrsmp(n, y, q, L)
%
% Bootstrap resampling for segmental bootstrapping of timeseries data
% n - number of resamples
% y - resampling matrix, l x q
% q - number of resamples in matrix
% L - the wanted size of the dataset

yl = size(y,1);
nsegs = ceil(L/yl);
outL = nsegs*yl;
X = zeros(outL, n);
for ii=1:n
   temp = zeros(outL,1);
   segs = randi(q, nsegs, 1); 
   for jj=1:nsegs
       temp(((jj-1)*yl+1):((jj)*yl)) = y(:,segs(jj));
   end
   X(:,ii) = temp;
end
