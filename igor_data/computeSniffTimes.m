function sweepStruct = computeSniffTimes(sweepStruct)
% function sweepStruct = computeSniffTimes(expFilename)
%
%
pb = 1; % to plot or not to plot
if (pb) figure; end
medFilt_len = 50; %10k sample rate, 5ms smoothing
hp_std = 300;
signal_name = 'Thermocouple'; %probably always Thermocouple or Pressure
thresh_sd = 4; %number of standard deviations to be considered a sniff
%filter parameters
sniff_filt_len = 100e-3; %gaussian to produce a sniff frequency vector, std in sec
hp_filt_len = 5e-3; %filter length for noise estimation, empirical, in s
peak_window = 60e-3; %window for detecting sniff peaks, necessary for peak finding algorithm but sets
                     % upper bound on detected frequency.

%sweepStruct = readIgorH5file(expFilename, signal_name, medFilt_len, hp_std);

for ii=1:length(sweepStruct)
    trace = sweepStruct(ii).value_filt;
    dt = sweepStruct(ii).dt;
    %make filter kernels
    %gaussian filtering to produce a sniff frequency vector
    %standard  deviation of the gaussian smoothing [in 1/10 ms]
    stdg = round(sniff_filt_len / dt); %2000 samples is a gaussian with std 200ms
    filtx = -3*stdg:1:3*stdg;
    filty = normpdf(filtx,0,stdg);
    filty = filty/sum(filty);
    % want to come up with a filter that gets rid of all of the noise, being a little overzealous.
    stdg_hp = round(hp_filt_len / dt); %this is a nice length, empirical, in ms
    filtx_hp = -3*stdg:1:3*stdg;
    filty_hp = normpdf(filtx_hp,0,stdg_hp);
    filty_hp = filty_hp/sum(filty_hp);
    filt_shift = floor(length(filtx_hp)/2);
    
    % need to figure out a way of getting a decent noise estimate in order to come up with a meaningful 
    % threshold for detecting peaks.
    filt_trace = conv2(trace, filty_hp(:), 'valid');
    filt_inds = (1:length(filt_trace)) + filt_shift;
    noise = trace(filt_inds) - filt_trace;
    noise_sd = std(noise(:)); % this is the value we wanted
    
    [sniff_delta, sniffs] = findPeaks(trace, noise_sd * thresh_sd, -1, peak_window / dt);
    %[sniff_delta, sniffs] = findPeaks(trace, abs(nanmin(trace))*.3, -1);
    sniff_bool = logical(sniff_delta);
    sniffy = trace(sniff_bool);
    sniffx = sweepStruct(ii).time(sniff_bool);
    if pb
        plot(sweepStruct(ii).time, trace, 'k-'); hold on;
        plot(sniffx, sniffy, 'ro');
    end
    sweepStruct(ii).sniffVect = sniff_bool;
    
    filtbuff = conv2(double(sniff_bool), filty(:), 'same');
    sniff_freq = filtbuff;
    %normalize to the number of sniffs in the trial as done for spikes in Gabbiani et al. 1999
    sniff_freq = sniff_freq /((sum(sniff_freq))*sweepStruct(ii).dt);
    sweepStruct(ii).sniffFreq = sniff_freq *sum(sweepStruct(ii).sniffVect);
    
    % Now, let's calculate the frequencies doing a reciprocal ISI strategy
    sniffi = find(sniff_bool);
    sniffISI = diff(sniffi).*sweepStruct(ii).dt;
    sniffInv = 1./sniffISI;
    isit = sweepStruct(ii).time(sniffi(1):sniffi(end));
    isir = zeros(length(sniffi(1):sniffi(end)),1);
    for jj=1:(length(sniffi)-1)
        start = sniffi(jj) - sniffi(1) +1;
        endi = sniffi(jj+1) - sniffi(1);
        isir(start:endi) = sniffInv(jj);
    end
    sweepStruct(ii).ISIrate_t = isit;
    sweepStruct(ii).ISIrate = isir;
end

    