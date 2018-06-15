% Logistic regression of the binary turning data

resampleSniffData;

%overall probability of a turn in the dataset
global_turnp = squeeze(sum(resamp_bturn,1)) ./ size(resamp_bturn,1);
global_mturnp = mean(global_turnp, 2);

% try some logistic regression of the variables
% First, just gonna plot the relevant vars across the range that we're considering
%for kk=1:rep
samp_dND = sniffData.dist_diffs(resamp_idx(:,3,1), :);
samp_ND = sniffData.sniffPos(resamp_idx(:,3,1), :);
bturn = resamp_bturn(:,3,1);
figure;
title ('Relationships of variables to turn prob if \Delta ND uniformly selected');
for i = 1:3
    subplot(2,4,i);
    plot(samp_dND(:,i), .5*bturn + .11,'k.'); hold on;
    [N,bin] = histc(samp_dND(:,i), edges);
    for j=1:length(N)-1
        tp(j) = sum(bturn(bin==j)) ./ N(j);
    end
    xval = (edges(1:end-1) + edges(2:end)) / 2;
    stairs(xval, tp, 'r');
    xlabel('\Delta ND'); ylabel('Turn Yes/No');
    ylim([0.1 0.62]);
    xlim([-9 9]);
end
%end 

edges = -15:15; %binning edges
for i = 1:4
    subplot(2,4,i+4); hold on;
    [N,bin] = histc(samp_ND(:,i), edges);
    for j=1:length(N)-1
        tp(j) = sum(bturn(bin==j)) ./ N(j);
    end
    xval = (edges(1:end-1) + edges(2:end)) / 2;
    stairs(xval, tp, 'r');
    plot(samp_ND(:,i), .5*bturn + .11,'k.');
    xlabel('ND'); ylabel('Turn Yes/No');
    ylim([0.1 0.62]);
    xlim([min(edges) max(edges)]);
end




% rsDND = NaN*zeros(50*sampN*(length(edges)-1),1);
% rsND = NaN*zeros(50*sampN*(length(edges)-1),1);
% for jj = 1:size(sniffData.dist_diffs,2)
%     for kk=1:rep
%         turnDND = [];
%         for ii = 1:(nbins-1)
%             [turnb, idx] = datasample(sniffData.bturn(bin == ii), sampN); %sample uniformly based on binning above
%             turnp(ii,jj,kk) = sum(turnb);
%             ti = (turnb == 1);
%             NDsamp =sniffData.sniffPos(bin == ii, jj+1); %corresponding NDs
%             NDsamp = NDsamp(idx);
%             dNDsamp = sniffData.dist_diffs(bin==ii,jj); 
%             dNDsamp = dNDsamp(idx);
%             turnDND = cat(1, turnDND, dNDsamp(ti));
%             if(kk ==1 && jj==3)
%                 toti = (kk)*sampN*(ii-1) + (1:sampN);
%                 rsDND(toti) = dNDsamp;
%                 rsND(toti) = NDsamp;
%             end
%         end
%         turnDNDcell{jj,kk} = turnDND;
%     end
% end
% turnp(isnan(turnp)) = 0; %replace possible NaNs
% mturnp = nanmean(turnp,3);
% rsDND = rsDND(~isnan(rsDND));
% rsND = rsND(~isnan(rsND));

% Figure 1 - Distance changes, turning probabilities
% ---------------------------------------------------
% figure;
% subplot(1,3,1);
% [N, bins] = histc(rsDND, edges);
% stairs(edges, N);
% hold on;
% [N_ND, bins] = histc(rsND, edges);
% stairs(edges, N_ND);
% xlabel('\Delta Distance from Trail (mm)');
% ylabel('Counts');
% xlim([-10 10]);
% 
% subplot(1,3,2);
% %bar(edges(1:end-1),turnp, 'histc');
% stairs(edges',mturnp./sampN);
% xlabel('\Delta Distance from Trail (mm)');
% ylabel('Turning Probability');
% xlim([-10 10]);
% ylim([0 .5]);