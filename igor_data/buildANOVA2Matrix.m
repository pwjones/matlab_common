
minelem = 0;
for kk = 1:numel(all_mags_dd)
    if kk==1
        minelem = numel(all_mags_dd{kk});
    else
        minelem = min(minelem, numel(all_mags_dd{kk}));
    end
end

testM = zeros(size(all_mags_dd,1)*minelem, size(all_mags_dd,2));
for ii = 1:size(all_mags_dd,1)
    for jj = 1:size(all_mags_dd,2)
        tempi = randperm(numel(all_mags_dd{ii,jj}), minelem);
        temp = all_mags_dd{ii,jj}(tempi);
        rowStart = ( (ii-1)*minelem ) + 1;
        rowEnd = (ii)*minelem;
        testM(rowStart:rowEnd,jj) = temp;
    end
end
    