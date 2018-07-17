%resampleSniffData
rep = 500;
sampN = 500;
nback = 3;
edges = -9:9; %binning edges
pos_edges = -15:15;
% resampling the sniffing data in order to create a uniform sampling across the variables that I want to look
% at turning probability relative to. First bin the variables and then randomly sample with replacement from each
% bin.

%initializing vars
tp=[];
turnp = zeros(length(edges), nback, rep);
ND_resamp = zeros(length(edges), nback, rep);
deltaND_resamp = NaN*zeros(length(edges), nback, rep);
turnDNDcell = cell(size(sniffData.dist_diffs,2), rep);
resamp_idx = []; resamp_bturn=[];
for kk=1:rep
    for jj = 1:size(sniffData.dist_diffs,2)
        %[N,bin] = histc(sniffData.sniffPos(:,jj+1), pos_edges);
        [N,bin] = histc(sniffData.dist_diffs(:,jj), edges);
        nbins = length(N);
        samp_idx = [];
        for ii = 1:(nbins-1)
            bin_idx = find(bin ==ii);
            [turnb, idx] = datasample(sniffData.bturn(bin == ii), sampN); %sample uniformly based on binning above
            samp_idx = cat(1, samp_idx, bin_idx(idx));
        end
        resamp_bturn(:,jj,kk) = sniffData.bturn(samp_idx); %turn booleans
        resamp_idx(:,jj,kk) = samp_idx; %the indices into the orginal dataset
    end
end
