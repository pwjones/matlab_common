%assume that sniffData is in the workspace
pd = fitdist(sniffData.followingSniffISI, 'gamma');
poissd = fitdist(sniffData.followingSniffISI, 'poisson');

[counts, binx] = hist(sniffData.followingSniffISI, 100);
countp = counts./numel(sniffData.followingSniffISI);
normcount = counts./max(counts);
minv = min(sniffData.followingSniffISI);
maxv = max(sniffData.followingSniffISI);
x = linspace(minv,maxv, 200);

gammay = pd.pdf(x);
gammay = gammay./max(gammay);
figure;
bar(binx, normcount, 'hist');
hold on;
plot(x, gammay);


%% Analysis to treat the sniff ISIs as a event vector
st = sniffData.followingSniffISI;
pt = findPeaks(st, .125, 1, 5);
figure; plot(st); hold on; plot(find(pt), st(logical(pt)), 'o');
ei = find(pt);
esep = diff(ei);
figure; hist(ei, 50);

1/mean(st)