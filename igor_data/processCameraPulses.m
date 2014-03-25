function camTrig = processCameraPulses(camTrig, vids, dbg)
% camTrig is the structure from the HF5 file.  Vids is the video strucure and
% can be empty if it's not ready. DBG is a debug plotting flag

% loop to find the frame indices - more explicitly the first sample of each high pulse in the trigger signal 
for ii=1:length(camTrig)
    camTrig(ii).value_bin = camTrig(ii).value > .5 * max(camTrig(ii).value); %Just threshold the trigger at half height
    ups = find(camTrig(ii).value_bin);
    jumps = find(diff(ups) > 1); %find breaks in the high pulses
    frames = [ups(1); ups(jumps+1)];
    camTrig(ii).frameInds = frames;
    camTrig(ii).nFrames = length(frames);
    if ii <= length(vids) % If the video analysis doesn't start at the begining of the video, then take into account
        camTrig(ii).vidStartFrame = vids(ii).frameRange(1);
        camTrig(ii).frameRange = vids(ii).frameRange(1):camTrig(ii).nFrames;
        camTrig(ii).nFrames = length(camTrig(ii).frameRange);
        camTrig(ii).analyzedFrameInds = camTrig(ii).frameInds(camTrig(ii).frameRange);
    else
        camTrig(ii).vidStartFrame = 1;
    end
    if dbg
        plot(camTrig(ii).time, camTrig(ii).value, 'k-', camTrig(ii).time(frames), camTrig(ii).value(frames), 'ro'); 
        hold on;
    end
end