function fitEventExponentials(miniTrace, dt)
% function fitEventExponentials(miniTrace, dt)
%
pb=1;
% Fits the exponential decay of the miniTrace
full_time = ((1:length(miniTrace)) - 1)*dt;
[miniPeak, peak_i] = nanmin(miniTrace);
% select the segment of the trace
slopeTrace = computeMiniSlope(miniTrace', dt);
decay = find(slopeTrace(peak_i:end) > 0);
decay_start = decay(1);
sec_deriv = diff(decay);
jumps = find(sec_deriv > 1);
if(~isempty(jumps))
   select_end = peak_i+jumps(1)+1; 
else
   select_end = peak_i+length(decay);
end
sel_trace = miniTrace(decay_start:select_end);
sel_time = full_time(decay_start:select_end);
if pb
    figure; hold on;
    plot(full_time, miniTrace);
    plot(sel_time, sel_trace, 'r');
end
t = (1:length(sel_trace))*dt;
exp_pred = @(a,xdata) a(1)*ones(size(t))+a(2)*exp(-t/a(3));

