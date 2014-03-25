distance_comp = cell(2,2); %want to do a reward/distractor, early/late comparison
for ii = 1:2
    trials = hist_trial_range(ii,:);
    for jj=1:length(trials)
        rdists = rew_dists{trials(jj)}; 
        nzi = rdists ~= 0; %these are single frame entries onto the trail
        rdists = rdists(nzi); %eliminate them
        
        ddists = distract_dists{trials(jj)};
        nzi = ddists ~= 0;
        ddists = ddists(nzi); 
        
        distance_comp{ii,1} = cat(1, distance_comp{ii,1}, rdists);
        distance_comp{ii,2} = cat(1, distance_comp{ii,2}, ddists);
    end
end
    
counts = cell(2,2); %want to do a reward/distractor, early/late comparison
xbins = 0:5:400;
[counts{1,1}] = hist(distance_comp{1,1}, xbins); counts{1,2} = hist(distance_comp{1,2},xbins);
counts{2,1} = hist(distance_comp{2,1},xbins); counts{2,2} = hist(distance_comp{2,2}, xbins);

% plotting things
xl = [-5 300];  %limits
yl = max([counts{1,1} counts{1,2}]); yl = [-(yl+5) yl+5];
dg = [0 .8 0]; %darker green
dr = [.8 0 0]; %darker red
figure;
% Early Trials
subplot(2,1,1); hold on;
bar(xbins, counts{1,1}, 'FaceColor', dg); hold on;
xmed = median(distance_comp{1,1});
line([xmed xmed], [0 yl(2)], 'Color', dg, 'LineStyle', '--'); 
text(xmed, yl(2)-10, ['median: ' num2str(xmed)], 'color', dg);

bar(xbins, -counts{1,2}, 'FaceColor', dr);
xmed = median(distance_comp{1,2});
line([xmed xmed], [0 yl(1)], 'Color', dr, 'LineStyle', '--'); 
text(xmed, yl(1)+10, ['median: ' num2str(xmed)], 'color', dr);
title('First 8 Trials');
set(gca, 'TickDir', 'out');
xlim(xl); ylim(yl);
% Late Trials
yl = max([counts{2,1} counts{2,2}]);  yl = [-(yl+5) yl+5];
subplot(2,1,2); hold on;
bar(xbins, counts{2,1}, 'FaceColor', dg); hold on;
xmed = median(distance_comp{2,1});
line([xmed xmed], [0 yl(2)], 'Color', dg, 'LineStyle', '--'); 
text(xmed, yl(2)-10, ['median: ' num2str(xmed)], 'color', dg);

bar(xbins, -counts{2,2}, 'FaceColor', dr); hold on;
xmed = median(distance_comp{2,2});
line([xmed xmed], [0 yl(1)], 'Color', dr, 'LineStyle', '--'); 
text(xmed, yl(1)+10, ['median: ' num2str(xmed)], 'color', dr);
xlim(xl); ylim(yl);
title('Last 8 Trials');
xlabel('Following Distance (px)', 'FontSize', 14);
ylabel('# Instances','FontSize', 14);
set(gca, 'TickDir', 'out');
