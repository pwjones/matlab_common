[meanSniffRate, sniffRate, respTrace, sniffs] = getSniffRates(exp, 2,300:1050, 0);
time = ((1:length(sniffRate))-1)*exp.resp(2).dt;
sniffi = find(sniffs);
rt = (respTrace - mean(respTrace));
rt = rt./(max(rt)-min(rt))*10 + 8;

figure; 
plot(time, rt, 'k'); 
hold on; 
plot(time, sniffRate, 'r', 'LineWidth', 2);
plot(time(sniffi), rt(sniffi), 'rx', 'MarkerSize', 12);
xlabel('Time (sec)');
ylabel('Respiration Frequency (Hz)');
xlim([1 5]);
ylim([3 15]);
set(gca, 'TickDir', 'out');