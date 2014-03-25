% plotRespirationTraces.m

if exist('exp', 'var')
    figure;
    for ii = 1:length(exp.resp)
        subplot (3, ceil(length(exp.resp)/3), ii);
        si = exp.resp(ii).sniffVect;
        plot(exp.resp(ii).time, exp.resp(ii).value_filt, 'k'); hold on;
        plot(exp.resp(ii).time(si), exp.resp(ii).value_filt(si), 'ro');  
        plot(exp.resp(ii).time, exp.resp(ii).sniffFreq * 10, 'b');
        
    end
    
else
    disp('The variable EXP doesnt exist in the workspace');
end