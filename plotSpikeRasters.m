function rh = plotSpikeRasters(axis_h, time, spikes, position)
% function rh = plotSpikeRasters(ah, time, spikes, position)
%
% Make a set of spike rasters from a delta function spike vector
 
ntrials = size(spikes,1);
spksize = size(spikes);
rh = zeros(nansum2(nansum2(spikes)), 1);
count = 1;
for i=1:ntrials
    y_vect = [(i) (i+0.9)] + position;
    spki = find(spikes(i,:) > 0.5); %spikes are 1, non spikes are 0/NaN
    if (~isempty(spki))
        for m = 1:length(spki)
          %spk_inds_m = spk_inds(m);
          t_val = time(spki(m));
          rh_temp = line('Parent',axis_h,'XData',[t_val t_val],...
              'YData',y_vect,'Tag','single_spike', 'LineWidth', 1);
          rh(count) =rh_temp;
          count = count + 1;
        end;
    end;
end
