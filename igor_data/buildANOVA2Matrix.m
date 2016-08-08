%find the minimum number of repetitions to use that number for making the 
%ANOVA matrix
minelem = 0;
for kk = 1:numel(all_mags_dd)
    if kk==1
        minelem = numel(all_mags_dd{kk});
    else
        minelem = min(minelem, numel(all_mags_dd{kk}));
    end
end
% build matrix of #Rows = DistanceBins * Reps, #Columns = Sniff Groups 
testM = zeros(size(all_mags_dd,1)*minelem, size(all_mags_dd,2));
for ii = 1:size(all_mags_dd,2)
    for jj = 1:size(all_mags_dd,1)
        tempi = randperm(numel(all_mags_dd{jj,ii}), minelem);
        temp = all_mags_dd{jj,ii}(tempi);
        rowStart = ( (jj-1)*minelem ) + 1;
        rowEnd = (jj)*minelem;
        testM(rowStart:rowEnd,ii) = temp;
    end
end
    