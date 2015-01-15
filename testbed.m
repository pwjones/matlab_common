
if ~isempty(ctl_trials{jj})
    epochi = find(trialCount >= ctl_trials{jj}(1) & trialCount <= ctl_trials{jj}(end));
    data = posTS.data(epochi);
end
