% dip testing the resampled data
sz = size(turnDNDcell);
dip = NaN*zeros(sz);
dip_pval = NaN*zeros(sz);
for ii=1:sz(1)
    parfor jj=1:sz(2)
        [dip(ii,jj),dip_pval(ii,jj)] = HartigansDipSignifTest(turnDNDcell{ii,jj}, 500);
    end
end
