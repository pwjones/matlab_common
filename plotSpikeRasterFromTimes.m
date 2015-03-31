function rh = plotSpikeRasterFromTimes(axis_h, times, position)
% function rh = plotSpikeRasters(ah, time, spikes, position)
%
% Make a set of spike rasters from a delta function spike vector
 
rh = zeros(length(times(:)), 1);
count = 1;

y_vect = [0 0.9] + position;
if (~isempty(times))
    for m = 1:length(times)
      rh_temp = line('Parent',axis_h,'XData',[times(m) times(m)],...
          'YData',y_vect,'Tag','single_spike', 'LineWidth', 1);
      %rh_temp = plot(axis_h, [t_val t_val], y_vect);
      rh(count) =rh_temp;
      count = count + 1;
    end;
end;
