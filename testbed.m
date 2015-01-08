%% Plotting the mean shifts and computing significance
figure;
pb = 0;
histx = -20:.5:20;
mean_turn_dists = squeeze(nanmean(turn_dists));
h = NaN*zeros(size(mean_turn_dists));
for ii = 2:size(turn_dists,3)
    for jj = 1:2
        for kk = 1:nMice
            x1 = turn_dists(:,jj,1,kk);
            x1 = x1(~isnan(x1));
            x2 = turn_dists(:,jj,ii,kk);
            x2 = x2(~isnan(x2));
            if pb
               figure;
               y1 = histc(x1, histx);
               y2 = histc(x2, histx);
               bar(histx', [y1, y2]);
            end
            
            if ~isempty(x1) && ~isempty(x2)
                [p,h] = ranksum(x1,x2);

                if (h == 1)
                    plot(ii, mean_turn_dists(jj,ii,kk) - mean_turn_dists(jj,1,kk), 'ko', 'MarkerFaceColor', 'k');
                else
                    plot(ii, mean_turn_dists(jj,ii,kk) - mean_turn_dists(jj,1,kk), 'ko');
                end
                hold on;
            end
        end
    end
end