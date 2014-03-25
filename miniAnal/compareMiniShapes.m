clear all;
% compareMiniShapes
ind_pb = 0;

dirname = '~pwjones/glab/data/lgmd_vc/090615/';
%cc_fname = 'man_8_mini.mat';
%vc_fname = 'man_10_mini.mat';
%dirname = '~pwjones/glab/data/lgmd_vc/090616/';
%cc_fname = 'man_5_mini.mat';
%vc_fname = 'man_6_mini.mat';
dirname = '~pwjones/glab/data/lgmd_vc/processed/minis/';
%cc_fname = '090508_1-92_mini.mat';
%vc_fname = '090508_1_1-78_mini.mat';
cc_fname = '090501_1_1-34_mini.mat';
%vc_fname = '090501_1-25_mini.mat';
%vc_fname = '090501_1-37_mini.mat';
%vc_fname = '071211_1-33_mini.mat';
cc_fname = '090508_1-92_mini.mat';
%vc_fname = '090508_1_1-78_mini.mat';
% load the files
cc = load([dirname cc_fname]);
cc = cc.handles;
vc = load([dirname vc_fname]);
vc = vc.handles;

% a little adaptation for mat files that have minis from multiple files
if (isfield(cc, 'prev') && length(cc.prev) > 0)
    [cc.minis, cc.miniSlopes, cc.miniHeights] = returnAllMinis(cc);
    cc.numMinis = sum([cc.prev.numMinis]) + cc.numMinis;
end
if (isfield(vc, 'prev') && length(vc.prev) > 0)
    [vc.minis, vc.miniSlopes, vc.miniHeights] = returnAllMinis(vc);
    vc.numMinis = sum([vc.prev.numMinis]) + vc.numMinis;
end

% exclude minis that somehow have positive heights
includei = find(vc.miniHeights < 0);
vc.miniHeights = vc.miniHeights(includei);
vc.minis = vc.minis(:,includei);
vc.miniSlopes = vc.miniSlopes(:,includei);
nexcluded = vc.numMinis - length(includei);
disp(sprintf('Excluded %i minis in the analysis', nexcluded));
vc.numMinis = length(includei);


mean_cc = nanmean2(cc.minis,2);
cc_len = length(mean_cc);
mean_cc = mean_cc - mean(mean_cc(1:floor(cc_len/2)));
min_cc = min(mean_cc);
norm_cc = mean_cc./-min_cc;
mean_vc = nanmean2(vc.minis,2);
mean_vc = mean_vc - mean(mean_vc(1:floor(cc_len/2)));
min_vc = min(mean_vc);
norm_vc = mean_vc./-min_vc;

% Do the CC (voltage trace "minis")
minis = zeros(cc.numMinis, cc.miniAnalWindow/cc.dt2) * NaN;
%windowVect = 1:(cc.miniAnalWindow/cc.dt2);
%windowVect = windowVect - floor((length(windowVect)/2));
windowLen = cc.miniAnalWindow/cc.dt2;
windowVect = (0:windowLen) - windowLen/2;
miniTime = ((1:length(windowVect))-1) * cc.dt2;
%handles = detectMiniOnsets(handles);

minis = cc.minis;
for i = 1:cc.numMinis
    neg_slopes = find(cc.miniSlopes(:,i) < 0);
    %minis(i, :) = cc.trace2_sub(windowVect+cc.miniPeaks(i));
    %minis(i, :) = cc.trace2(windowVect+cc.miniOnsets(i)) - cc.trace2(cc.miniOnsets(i));
    
    % onset is aligned at half. Going to average till the EPSC upslope ends, then it is NaN.  
    slopeN = size(cc.miniSlopes,1);
    last_half = floor(slopeN/2)+1:slopeN; 
    neg_slopes = intersect(last_half, find(cc.miniSlopes(:,i) < 0));
    pos_slopes = setdiff(last_half, neg_slopes); %these are the upslopes
    jumps = find(diff(pos_slopes) > 1); %non-contiguous segments
    if isempty(jumps)
        change_ind = slopeN;
    else
        change_ind = pos_slopes(jumps(1))+1;
    end
    if ind_pb
        figure; hold on;
        plot(miniTime, minis(:,i), 'k');
        plot(miniTime(neg_slopes), minis(neg_slopes, i), 'b');
        plot(miniTime(pos_slopes), minis(pos_slopes, i), 'g');
        plot(miniTime(change_ind), minis(change_ind, i), 'gx');
        plot(miniTime, cc_meanMini, 'r', 'Linewidth', 2);
    end
    minis(change_ind:end, i) = NaN;
end
cc_meanMini = nanmean2(minis,2);


minis = zeros(size(vc.minis,1), vc.numMinis) * NaN;
%windowVect = 1:(vc.miniAnalWindow/vc.dt2);
%windowVect = windowVect - floor((length(windowVect)/2));
%miniTime = ((1:length(windowVect))-1) * vc.dt2;
%handles = detectMiniOnsets(handles);

for i = 1:vc.numMinis
    neg_slopes = find(vc.miniSlopes(:,i) < 0);
    %minis(i, :) = vc.trace2_sub(windowVect+vc.miniPeaks(i));
    %minis(i, :) = vc.trace2(windowVect+vc.miniOnsets(i)) - vc.trace2(vc.miniOnsets(i));
    minis(:,i) = vc.minis(:,i);
    
    % onset is aligned at half. Going to average till the EPSC upslope ends, then it is NaN.  
    slopeN = size(vc.miniSlopes,1);
    last_half = floor(slopeN/2)+1:slopeN; 
    neg_slopes = intersect(last_half, find(vc.miniSlopes(:,i) < 0));
    pos_slopes = setdiff(last_half, neg_slopes); %these are the upslopes
    jumps = find(diff(pos_slopes) > 1); %non-contiguous segments
    if isempty(jumps)
        change_ind = slopeN;
    else
        change_ind = pos_slopes(jumps(1))+1;
    end
    if ind_pb
        figure; hold on;
        plot(miniTime, minis(:,i), 'k');
        plot(miniTime(neg_slopes), minis(neg_slopes, i), 'b');
        plot(miniTime(pos_slopes), minis(pos_slopes, i), 'g');
        plot(miniTime(change_ind), minis(change_ind, i), 'gx');
        plot(miniTime, vc_meanMini, 'r', 'Linewidth', 2);
    end
    minis(change_ind:end, i) = NaN;
end
vc_meanMini = nanmean2(minis,2);



mean_cc = cc_meanMini;
cc_len = length(mean_cc);
mean_cc = mean_cc - mean_cc(floor(cc_len/2));
min_cc = min(mean_cc);
norm_cc = mean_cc./-min_cc;

mean_vc = vc_meanMini;
mean_vc = mean_vc - mean_vc(floor(cc_len/2));
min_vc = min(mean_vc);
norm_vc = mean_vc./-min_vc;
norm_cc = mean_cc ./ min_cc * min_vc; %Actually normalize the CC trace to 1, the multiply by the VC height 

miniTime = miniTime - miniTime(floor(cc_len/2));
figure; hold on;
ph(1) = plot(miniTime, mean_vc, 'k', 'LineWidth',2);
ph(2) = plot(miniTime, norm_cc, 'Color', [.6 .6 .6], 'LineWidth',2, 'Linestyle', '--');
legend({'EPSC','Inverted/Scaled EPSP'});
title(sprintf('Mean event shape for experiment %s', dirname));
xlabel('Time (ms)', 'FontSize', 14);
ylabel('Spont EPSC amplitude (nA)', 'FontSize', 14);
plot([miniTime(1), miniTime(end)], [0 0], 'k');
xlim([-1 6]);

%let's figure out the integral of the mean VC trace
vc_integral = quad(@evalMini, 0, 6, [], [], miniTime(:), vc_meanMini(:))
cc_integral = quad(@evalMini, 0, 6, [], [], miniTime(:), cc_meanMini(:))

% lets fit the mini traces
%fitEventExponentials(vc_meanMini(:), vc.dt2);
