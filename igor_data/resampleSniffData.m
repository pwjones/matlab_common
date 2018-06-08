%resampleSniffData
rep = 500;
sampN = 500;
nback = 3;
edges = -9:9; %binning edges

% Figure 1 - Distance changes, turning probabilities
% ---------------------------------------------------
%initializing vars
turnp = zeros(length(edges), nback, rep);
ND_resamp = zeros(length(edges), nback, rep);
deltaND_resamp = NaN*zeros(length(edges), nback, rep);
rsDND = NaN*zeros(50*sampN*(length(edges)-1),1);
rsND = NaN*zeros(50*sampN*(length(edges)-1),1);
turnDNDcell = cell(size(sniffData.dist_diffs,2), rep);
for jj = 1:size(sniffData.dist_diffs,2)
    %[N,bin] = histc(sniffData.sniffPos(:,jj+1), edges);
    [N,bin] = histc(sniffData.dist_diffs(:,jj), edges);
    nbins = length(N);
    for kk=1:rep
        turnDND = [];
        for ii = 1:(nbins-1)
            [turnb, idx] = datasample(sniffData.bturn(bin == ii), sampN); %sample uniformly based on binning above
            turnp(ii,jj,kk) = sum(turnb);
            ti = (turnb == 1);
            NDsamp =sniffData.sniffPos(bin == ii, jj+1); %corresponding NDs
            NDsamp = NDsamp(idx);
            dNDsamp = sniffData.dist_diffs(bin==ii,jj); 
            dNDsamp = dNDsamp(idx);
            turnDND = cat(1, turnDND, dNDsamp(ti));
            if(kk ==1 && jj==3)
                toti = (kk)*sampN*(ii-1) + (1:sampN);
                rsDND(toti) = dNDsamp;
                rsND(toti) = NDsamp;
            end
        end
        turnDNDcell{jj,kk} = turnDND;
    end
end
turnp(isnan(turnp)) = 0; %replace possible NaNs
mturnp = nanmean(turnp,3);
rsDND = rsDND(~isnan(rsDND));
rsND = rsND(~isnan(rsND));

figure;
subplot(1,3,1);
[N, bins] = histc(rsDND, edges);
stairs(edges, N);
hold on;
[N_ND, bins] = histc(rsND, edges);
stairs(edges, N_ND);
xlabel('\Delta Distance from Trail (mm)');
ylabel('Counts');
xlim([-10 10]);

subplot(1,3,2);
%bar(edges(1:end-1),turnp, 'histc');
stairs(edges',mturnp./sampN);
xlabel('\Delta Distance from Trail (mm)');
ylabel('Turning Probability');
xlim([-10 10]);
ylim([0 .5]);