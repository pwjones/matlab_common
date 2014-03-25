for ii = 1:length(sweeps)
    filt_trace{ii} = medfilt1(double(sweeps(ii).value), 20);
end

%%

figure; hold on;
for ii = 1:length(sweeps)
    plot(sweeps(ii).time, sweeps(ii).value);
end

%%

figure; hold on;
for ii = 1:length(sweepStruct)
    plot(sweepStruct(ii).time, sweepStruct(ii).sniffFreq);
end

%%

figure; hold on;
for ii = 1:length(sweepStruct)
    plot(sweepStruct(ii).ISIrate_t, sweepStruct(ii).ISIrate);
    plot(sweepStruct(ii).time, sweepStruct(ii).sniffFreq, 'r');
end

%%

figure; hold on;
for ii = 1:length(sweeps)
    plot(pulses(ii).time, pulses(ii).value);
end

%%
sniff_frames = zeros(length(sniff_times),1);
for ii=1:length(sniff_times)
    sniff_frames(ii) = find(frame_t >= sniff_times(ii), 1, 'first');
end

%%
% Let's get the frames when the animal sniffed
for ii=1:length(exp.resp)
    frame_t = exp.camTrig(ii).time(exp.camTrig(ii).frameInds);
    sniff_times = exp.resp(ii).time(exp.resp(ii).sniffVect);
    sniff_frames = zeros(length(sniff_times),1);
    for jj=1:length(sniff_times)
        sniff_frames(jj) = find(frame_t >= sniff_times(jj), 1, 'first');
    end
    exp.resp(ii).sniff_times = sniff_times;
    exp.resp(ii).sniff_frames = sniff_frames;
    exp.resp(ii).sniff_pos = exp.vids(ii).nosePos(sniff_frames,:);
end