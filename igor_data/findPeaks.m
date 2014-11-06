function [peakTimes, eventVects] = findPeaks(signal, thresh, sign, peak_window)
% function [peakTimes, eventVects] = findPeaks(signal, thresh, sign, peak_window)
% 
% signal - the signal you want to detect peaks in
% thresh - the threshold value for peak detection
% sign - the direction, 1 for positive, -1 for negative, of peaks
% peak_window - the number of samples smaller than which you can't have multiple peaks.

if sign == -1
    signal = -signal;
end
eventVects = zeros(1, 2*peak_window+1);

thresdw = thresh/2;

n_pts = length(signal);
peakTimes = zeros(size(signal));
event_x = -peak_window:peak_window;

done = 0;
spknum = 0;
i = 1;
while ( ~done )
    while ( (i<n_pts) & (signal(i)<thresh) )
        i = i + 1;
    end;
    
    if ( i >= n_pts )
        %end of trace
        done = 1;
    else
        %found a threshold crossing; now find the peak value within
        % the given window
        spike_start = i;
        peak_val = signal(i);
%         while ( (i<(spike_start+peak_window)) & (i <= n_pts)  ...
%                 &  (signal(i)>thresdw) )
%             if ( signal(i) > peak_val )
%                 peak_val = signal(i);
%                 spike_peak = i;
%             end;
%             i = i + 1;
%         end;
        win_end = min(n_pts, i+peak_window); %look ahead
        [peak_val, spike_peak] = max(signal(i:win_end)); 
        spike_peak = spike_peak+i;
        peakTimes(spike_peak) = 1;
        spknum = spknum+1;
        window_min = max(spike_peak - peak_window, 1);
        window_max = min(spike_peak + peak_window, n_pts);
        win_len = length(window_min:window_max);%window_max-window_min+1;
        %center = peak_window+1;
        %eventVects(spknum, ((peak_window-floor(win_len/2)):(peak_window+floor(win_len/2)))+1) = signal(window_min:window_max);
%         if ( (i <= spike_start+peak_window) & (i <= n_pts) & ...
%                 (signal<= thresdw) )
%             %really found a spike
%             %disp('adding a spike');
%             peakTimes(spike_peak) = 1;
%             spknum = spknum+1;
%             eventVects(spknum,:) = signal(spike_peak-peak_window:spike_peak+peak_window);
%         end;
        %i = spike_start+peak_window; %this allows events to be spaced too closely
        i = spike_peak+(peak_window/2);
    end;
end;


