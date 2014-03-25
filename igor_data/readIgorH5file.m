function data = readIgorH5file(filename, signalName, medfilt_len, hp_std)
% function data = readIgorH5file(filename, signalName)
% 
% Reads an HDF5 exported Igor experiment file for the traces (sweeps) of a
% particular signal name (signalName).  Assumes structure as given by
% exporting an experiment file (*.pxp) from Igor 6 using Rick Gerkin's
% RecordingArtist package for experimental control
%
% Both 

basepath = ['/root/' signalName];
if (~exist(filename, 'file'))
    error('File given does not exist');
end
infoT = h5info(filename, basepath);

sweepNames = cell(length(infoT.Datasets), 1);
for ii=1:length(infoT.Datasets)
    sweepNames{ii} = infoT.Datasets(ii).Name;
end
isSweep = strncmp('sweep', sweepNames, 5); %we are really only concerning ourselves with data labeled 'sweep'

jj=1;
for ii=1:length(infoT.Datasets)
    temp = [];
    temp.name = infoT.Datasets(ii).Name;
    if isSweep(ii)
        temp.value = double(h5read(filename, [basepath '/' temp.name]));
        temp.nSamp = length(temp.value);
        timeField = findh5FieldNumber(infoT.Datasets(ii).Attributes, 'IGORWaveNote');
        timeStr = infoT.Datasets(ii).Attributes(timeField).Value; %in the format TIME=xxx.xxx
        temp.ts = parseTimeStampStr(timeStr);
        if jj==1 firstTS = temp.ts; end
        % timestamps are really big if you don't subtract first time - also
        % adjusting to seconds from microseconds (timestamp is microseconds since boot)
        temp.ts = (temp.ts - firstTS) * 1e-6; 
        waveField = findh5FieldNumber(infoT.Datasets(ii).Attributes, 'IGORWaveScaling');
        waveScaling = infoT.Datasets(ii).Attributes(waveField).Value; 
        temp.dt = waveScaling(1,2); %value of sampling interval - in seconds
        
        temp.time = (((0:(length(temp.value)-1))*temp.dt)+temp.ts)'; %time vector for sweep relative to exp start
        data(jj) = temp;
        jj = jj+1;
    end
end

% because the string ordering throws off the actual trace order, need to
% reorder by timestamp
[~, oi] = sort([data.ts]);
data = data(oi);

% Now we should concatentate those sweeps that are actually continuously
% collected.
ts = [data.ts];
predTS = ts + [data.dt].*[data.nSamp];
pred_diff = ts(2:end) - predTS(1:end-1);
%cum_pred_diff = 
new_sweeps = [1 find(abs(pred_diff) > .5)+1]; %separate sweeps more than .5 second off
if length(new_sweeps) ~= length(data)
    for ii = 1:length(new_sweeps)
        first = new_sweeps(ii);
        if ii < length(new_sweeps)
            last = new_sweeps(ii+1) - 1;
        else
            last = length(data);
        end
        temp = [];
        temp.ts = data(first).ts;
        temp.dt = data(first).dt;
        temp.nSamp = sum([data(first:last).nSamp]);
        tv = [data(first:last).time];
        temp.time = tv(:);
        tv = [data(first:last).value];
        temp.value = tv(:);
        
        new_data(ii) = temp;
    end
    data = new_data;
end

% Now we can filter if we want, now that things are properly concatenated
if (medfilt_len || hp_std)
    data(1).value_filt = [];
    for ii =1:length(data)
        mf_samp_len = medfilt_len / data(ii).dt; %set the lengths in # of samples
        hp_samp_std = hp_std / data(ii).dt;
        if medfilt_len
            data(ii).value_filt = medfilt1(data(ii).value, mf_samp_len);
        end
        if hp_std
            baseline = gaussianFilter(data(ii).value_filt, hp_samp_std); 
            data(ii).value_filt = data(ii).value_filt - baseline;
        end
    end
end
