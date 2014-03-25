function varargout = miniAnal(varargin)
% MINIANAL M-file for miniAnal.fig
%      MINIANAL, by itself, creates a new MINIANAL or raises the existing
%      singleton*.
%
%      H = MINIANAL returns the handle to a new MINIANAL or the handle to
%      the existing singleton*.
%
%      MINIANAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MINIANAL.M with the given input arguments.
%
%      MINIANAL('Property','Value',...) creates a new MINIANAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before miniAnal_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to miniAnal_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help miniAnal

% Last Modified by GUIDE v2.5 07-Jul-2009 16:03:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @miniAnal_OpeningFcn, ...
                   'gui_OutputFcn',  @miniAnal_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before miniAnal is made visible.
function miniAnal_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to miniAnal (see VARARGIN)

% Choose default command line output for miniAnal
handles.output = hObject;

% Make some initial definitions and initialize some important vars
handles.dt1 = .05; %time between samples in msec
handles.dt2 = .05;
handles.intra_amp = 100;
handles.dv = 10/(2^16)*1000/handles.intra_amp; % (Volts/bit) * (mV/V) / amplification
handles.di = 10/(2^16)*1; % (Volts/bit) * 1nA/Volt sensitivity;
handles.medFilterWidth = 8; %number of samples to include in the median filter
%filter kernel for the low pass filtering to get a running baseline 
handles.stdg = 300; %300 samples is a gaussian with std 15ms
filtx = -3*handles.stdg:1:3*handles.stdg;
filty = my_normpdf(filtx,0,handles.stdg);
handles.filty = filty/sum(filty);
handles.riseTime = .5; %ms, event must rise in less time than this
handles.decayTime = 2;
handles.threshold = -.2; %nA, event must be more negative than this
handles.riseThreshold = -.1; %min nA rise in 2 ms
handles.slopeThreshold = -1; %min nA/sec rise
handles.onsetSlope = -.2;  %nA/sec rise
handles.miniAnalWindow = 20; %ms, length of the analysis window
handles.miniSearchWindow = 5; %ms, length of search window
handles.trace2_amp = 1; %amplification multiplier
handles.miniPeaks = [];
handles.numMinis = 0;
handles.numMinis_currFile = 0;
handles.minis = []; %while the 
handles.miniSlopes = [];
handles.miniHeights = [];
handles.miniPeaks = [];
handles.miniOnsets = [];
% handles to different plot items
handles.vmplot = [];
handles.implot = [];
handles.implot_filt = [];
handles.miniplot = [];
handles.thresholdline = [];
handles.rootdir = [];
handles.expdir = [];
handles.fname = [];
handles.currn = [];
handles.xlim = [0 500];
handles.traceType = 1;  %value for the traceType_menu callback, selects current(0) or voltage(1)
handles.invertTrace2 = 0;
handles.invertTrace1 = 0;
handles.format = 'int16';
handles.keepMinis = 0;
handles.prev = struct([]); % the structure in which to put previous files info, filenames, minilist, etc
% need to divide really long traces into segments and only plot a section
% in order to speed up the gui
handles.plotSection = 1;
handles.MAX_PLOT_LENGTH = 1e5;
handles.MAX_TRACE_TIME = 10000; %10 sec

%set up the limit windows for the display
set(handles.tmin_edit, 'String', num2str(handles.xlim(1)));
set(handles.tmax_edit, 'String', num2str(handles.xlim(2)));
set(handles.riseTime_edit, 'String', num2str(handles.riseTime));
set(handles.decayTime_edit, 'String', num2str(handles.decayTime));
set(handles.riseThreshold_edit, 'String', num2str(handles.riseThreshold));
set(handles.slopeThreshold_edit, 'String', num2str(handles.slopeThreshold));
set(handles.onsetSlope_edit, 'String', num2str(handles.onsetSlope));
set(handles.trace2amp_edit, 'String', num2str(handles.trace2_amp));
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes miniAnal wait for user response (see UIRESUME)
% uiwait(handles.miniAnal_fig);

% --- Outputs from this function are returned to the command line.
function varargout = miniAnal_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prompt  = {'Enter root directory:','Enter experiment directory:', ...
	   'Enter file base name:', 'Enter start protocol number:'};
title   = 'Initialization of thresholder';
lines= 1;

if ( ~isempty(handles.rootdir) )
  def_rootdir = handles.rootdir;
else
  def_rootdir = '~pwjones/glab/data/lgmd_vc/processed';
end;

if ( ~isempty(handles.expdir) )
  def_expdir = handles.expdir;
else 
  def_expdir = '090501';
end;

if ( ~isempty(handles.fname) )
  def_fname = handles.fname;
else
  def_fname = 'spont';
end;

if ( ~isempty(handles.currn) )
  def_filenum = num2str(handles.currn);
else
  def_filenum = '1';
end;
def     = {def_rootdir,def_expdir,def_fname,def_filenum};
answer  = inputdlg(prompt,title,lines,def);

if ( isempty(answer) )
  return;
end;

handles.rootdir = answer{1};
handles.expdir = answer{2};
handles.fname = answer{3};
handles.currn = str2num(answer{4});
guidata(hObject, handles);
handles.keepMinis = 0;
loadTrial(hObject, eventdata, handles);


% -------------------------------------------------------------
 function loadTrial(hObject, eventdata, handles)
% Loads the trial specified in the handles.structure

%First, save any data from previousl loaded files, if necessary
%The caller sets handles.keepMinis to decide this

% Now check if this is a file advancement and if there are already detected
% minis
% if(handles.numMinis ~= 0 & handles.keepMinis)
%     resp = questdlg('Do you want to keep the minis in the current file?');
%     if (strcmp('yes', resp))
%         handles.keepMinis = 1;
%     else
%         handles.keepMinis = 1;
%     end
% end
if (handles.keepMinis && handles.numMinis > 0)
    nprevfiles = length(handles.prev);
    prev.filename = handles.filename;
    prev.numMinis = handles.numMinis;
    prev.miniSlopes = handles.miniSlopes;
    prev.miniHeights = handles.miniHeights;
    prev.miniPeaks = handles.miniPeaks;
    prev.miniOnsets = handles.miniOnsets;
    prev.minis = handles.minis;
    prev.miniSlopes = handles.miniSlopes;
    prev.trace2_filt = handles.trace2_filt;
    prev.trace2 = handles.trace2;
    prev.time2 = handles.time2;
    prev.dt2 = handles.dt2;
    prev.trace1_filt = handles.trace1_filt;
    prev.trace1 = handles.trace1;
    prev.time1 = handles.time1;
    prev.dt1 = handles.dt1;
    
    if nprevfiles == 0 %to cover the case where this is the first prev file
        handles.prev = prev;
    else
        handles.prev(nprevfiles+1) = prev;
    end
    
    handles.numMinis = sum([handles.prev.numMinis]) + handles.numMinis;
    set(handles.prevMinis_text, 'string', num2str(sum([handles.prev.numMinis])));
    handles.numMinis = 0;
    set(handles.miniNum_text, 'string', num2str(handles.numMinis));
else
    handles.numMinis = 0;
    set(handles.miniNum_text, 'string', num2str(handles.numMinis));
end

% clear any previous data
handles.trace1 = [];
handles.trace1_filt = [];
handles.trace2 = [];
handles.trace2_filt = [];
handles.trace2_lp = [];
handles.time1 = [];
handles.time2 = [];
handles.time2_lp = [];
handles.miniPeaks = [];
handles.miniOnsets = [];
handles.miniHeights = [];
handles.miniSlopes = [];

%set vars to enable switching between analyzing vm and im data
imchan = 3;
vmchan = 2;
if handles.traceType == 1
    t1chan = vmchan;
    t2chan = imchan;
    t1calib = handles.dv;
    t2calib = handles.di;
else 
    t1chan = imchan;
    t2chan = vmchan;
    t1calib = handles.di;
    t2calib = handles.dv;
end

%load the first trace file
chan1file = sprintf('%s/%s/%s_%i_chan_%i',handles.rootdir,handles.expdir, ...
		    handles.fname, handles.currn, t1chan);
if ( exist(chan1file) ~= 2 )
    err_str = sprintf('unable to find %s',chan1file);
    errordlg(err_str,'Load trace error');
    return;
end;

fid = fopen(chan1file, 'r','l');
if ( fid == -1 )
    err_str = sprintf('error opening %s',chan1file);
    errordlg(err_str,'Load trace error');
    return;
end;
handles.trace1 = fread(fid,inf,handles.format);
max_samples = floor(handles.MAX_TRACE_TIME / handles.dt1);
if length(handles.trace1) > max_samples
    handles.trace1 = handles.trace1(1:max_samples);
end
fclose(fid);
if strcmp(handles.format, 'int16')
    handles.trace1 = handles.trace1 * t1calib; %convert A/D values to mV or nA
end
if (handles.invertTrace1) handles.trace1 = -handles.trace1; end
handles.trace1_filt = medfilt2(handles.trace1, [handles.medFilterWidth, 1]); %filter the trace2

%load the 2nd trace file - this is the one to be analyzed
chan2file = sprintf('%s/%s/%s_%i_chan_%i',handles.rootdir,handles.expdir, ...
        handles.fname, handles.currn, t2chan);

if ( exist(chan2file) ~= 2 )
    err_str = sprintf('unable to find %s',chan2file);
    errordlg(err_str,'Load intracellular trace error');
    return;
end;
fid = fopen(chan2file, 'r','l');
if ( fid == -1 )
    err_str = sprintf('error opening %s',chan2file);
    errordlg(err_str,'Load trace error');
    return;
end;
handles.trace2 = fread(fid, inf, handles.format);
max_samples = floor(handles.MAX_TRACE_TIME / handles.dt2);
if length(handles.trace2) > max_samples
    handles.trace2 = handles.trace2(1:max_samples);
end
fclose(fid);
handles.filename = chan2file; % update the current filename here, after potentially saving old info
if strcmp(handles.format, 'int16')
    handles.trace2 = handles.trace2 * t2calib; %convert A/D values to nA or mV
end
if (handles.invertTrace2) handles.trace2 = -handles.trace2; end
handles.trace2 = handles.trace2 * handles.trace2_amp; %additional amplification factor for user in GUI
if (strcmp(handles.format, 'double') && handles.traceType == 1) %These are preprocessed Im traces
    handles.trace2_filt = handles.trace2;
else
    handles.trace2_filt = medfilt2(handles.trace2, [handles.medFilterWidth, 1]); %filter the trace2
end
%come up with the proper time vectors
handles.time1 = (0:(length(handles.trace1)-1)) * handles.dt1;
handles.time2 = (0:(length(handles.trace2)-1)) * handles.dt2;
% The running avg time vector must be shorter by the half width of the
% filter
%handles.time2_lp = handles.time2(1:(length(handles.time2)-3*handles.stdg));
%gaussian smoothing
%filtbuff = fftfilt2(handles.filty,handles.trace2);
%correct for the time shift generated by the smoothing operation
%handles.trace2_lp = filtbuff((1+3*handles.stdg):end);

%redoing the longpass filtering step with filtfilt
handles.time2_lp = handles.time2;
handles.trace2_lp = filtfilt(handles.filty, 1, handles.trace2);
% also generate a running baseline subtracted trace
handles.trace2_sub = handles.trace2_filt(1:length(handles.time2_lp)) - handles.trace2_lp;

if (handles.time1(end) > handles.MAX_TRACE_TIME)
    over1 = find(handles.time1 > handles.MAX_TRACE_TIME);
    over2 = find(handles.time2 > handles.MAX_TRACE_TIME);
    over2_lp = find(handles.time2_lp > handles.MAX_TRACE_TIME);
    over1 = over1(1) - 1;
    over2 = over2(1) - 1;
    over2_lp = over2_lp(1) -1;
    handles.time1 = handles.time1(1:over1);
    handles.time2 = handles.time2(1:over2);
    handles.time2_lp = handles.time2_lp(1:over2_lp);
    handles.trace1 = handles.trace1(1:over1);
    handles.trace2 = handles.trace2(1:over2);
    handles.trace2_lp = handles.trace2_lp(1:over2_lp);
    handles.trace2_sub = handles.trace2_sub(1:over2_lp);
    
end

% plot the traces
if ( ~isempty(handles.vmplot) )
    set(handles.vmplot, 'XData', handles.time1, ...
        'YData', handles.trace1);
else
    handles.vmplot = line('Parent', handles.vm_axis, ...
                             'XData', handles.time1, ...
                             'YData', handles.trace1_filt, ...
                             'Color', 'b');
end;                   
set(handles.vm_axis, 'Xlim', handles.xlim, 'Ylim', [-Inf Inf]);

if ( ~isempty(handles.implot_filt) )
%     set(handles.implot, 'XData', handles.time2, ...
%         'YData', handles.trace2);
    set(handles.implot_filt, 'XData', handles.time2_lp, ...
        'YData', handles.trace2_sub);
    set(handles.implot_lp, 'XData', handles.time2_lp, ...
        'YData', handles.trace2_lp);
    set(handles.thresholdline, 'XData', [0 handles.MAX_TRACE_TIME], ...
        'YData', [handles.threshold handles.threshold]);
else
%     handles.implot = line('Parent', handles.im_axis, ...
%                              'XData', handles.time2, ...
%                              'YData', handles.trace2, ...
%                              'Color', 'k');
    handles.implot_filt = line('Parent', handles.im_axis, ...
                             'XData', handles.time2_lp, ...
                             'YData', handles.trace2_sub, ...
                             'Color', 'k', ...
                             'ButtonDownFcn', @addClickedEvent);
    handles.implot_lp = line('Parent', handles.im_axis, ...
                             'XData', handles.time2_lp, ...
                             'YData', handles.trace2_lp, ...
                             'Color', [.5 .5 .5]);
     handles.thresholdline = line('Parent', handles.im_axis, ...
                              'XData', [0 handles.MAX_TRACE_TIME], ...
                             'YData', [handles.threshold handles.threshold], ...
                             'Color', 'r', 'ButtonDownFcn', @start_move);
end; 
if ( ~isempty(handles.miniplot))
    set(handles.miniplot, 'xdata', [], 'ydata', []);
end
set(handles.fileNum_text, 'String', num2str(handles.currn));
set(handles.im_axis, 'Xlim',handles.xlim, 'Ylim', [-Inf Inf]);
autoSetYLimits(hObject, eventdata, handles);
handles = updateMiniTraces(hObject, handles);
handles = markMinis(gcbf, [], handles);
guidata(hObject,handles); %saves the altered handles data structure for fig h.


% --------------------------------------------------------------------
function handles = detectMinis(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Goes through the trace2 and detects the miniature events 
% It does this by defining a threshold deflection amount that the trace
% needs to change within a set amount of time.  The threshold is
% handles.threshold and the amount of time is handles.riseTime. 

handles.miniPeaks = [];
handles.miniOnsets = [];
handles.numMinis = 0;
correctedTrace = handles.trace2_sub;
riseSamples = handles.riseTime / handles.dt2;
windowLen = handles.miniSearchWindow/handles.dt2;
%windowVect = (0:windowLen) - windowLen/2;
windowVect = (0:windowLen);
analLen = handles.miniAnalWindow/handles.dt2;
analVect = (0:analLen) - analLen/2; 
% The dumb loop way of detecting threshold crossings, just stepping through
% the trace point by point, with a sliding window centered on the 
i = abs(analVect(1)) +1;
%i=1;
while i <= (length(correctedTrace) - windowLen/2)
    if(correctedTrace(i) <= handles.threshold)
        [localMin, localMinIndex] = min(correctedTrace(i+windowVect));
        localMinIndex = windowVect(localMinIndex)+i;
        startIndex = localMinIndex-riseSamples;
        if startIndex < 1
            startIndex = 1;
        end
        startVal = correctedTrace(startIndex);
        difft = correctedTrace(localMinIndex) - startVal;
        riseTimeMin = min(correctedTrace(startIndex:localMinIndex-1)); %within a rise time, it is a downward trace
        %found an event, checking for a min deflection and that there are
        %no intervening peaks
        [onsetInd, slopeTrace] = detectSingleMiniOnset(handles.trace2_filt(i+analVect), handles.dt2, handles.slopeThreshold);
        if (difft <= handles.riseThreshold && riseTimeMin > localMin && ~isnan(onsetInd)) 
            handles.miniPeaks(handles.numMinis +1) = localMinIndex;
            [onsetInd, slopeTrace] = detectSingleMiniOnset(handles.trace2_filt(i+analVect), handles.dt2, handles.onsetSlope);
            handles.miniOnsets(handles.numMinis+1) =  analVect(onsetInd)+i;
            handles.numMinis = handles.numMinis + 1;
            %i = localMinIndex + windowLen/2; %skip to the peak
            i = localMinIndex;
        end
    end
    %i = i+windowLen/2;
    i = i+1;
end
%uniquify the detected minis
%[b, m, n] = unique(handles.miniPeaks, 'first');
%handles.miniPeaks = b;
%handles.numMinis = length(handles.miniPeaks);
handles = updateMiniTraces(hObject, handles);
guidata(hObject, handles);

% --------------------------------------------------------------------
function handles = detectMinis2(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Goes through the trace2 and detects the miniature events 
% It does this by defining a threshold deflection amount that the trace
% needs to change within a set amount of time.  The threshold is
% handles.threshold and the amount of time is handles.riseTime. 

handles.miniPeaks = [];
handles.miniOnsets = [];
handles.numMinis = 0;
correctedTrace = handles.trace2_sub;
riseSamples = handles.riseTime / handles.dt2;
windowLen = handles.miniSearchWindow/handles.dt2;
%windowVect = (0:windowLen) - windowLen/2;
windowVect = (0:windowLen);
analLen = floor(handles.miniAnalWindow/handles.dt2);
analCenter = ceil(analLen/2);
analVect = (1:analLen) - analCenter;

handles.trace2slope = computeMiniSlope(handles.trace2_filt, handles.dt2, 25);
posi = find(handles.trace2slope >= 0);
negi = find(handles.trace2slope < 0);
neg_jump = find(diff(negi) > 1); %find breaks in positive slopes that are zero crossings - local maxima
minima = negi(neg_jump);
pos_jump = find(diff(posi) > 1); %find breaks in positive slopes that are zero crossings - local maxima
maxima = posi(pos_jump);
max_pairs = NaN*zeros(size(minima));
for i=1:length(minima) %now find the local maxima that precede the local minima
    max_diff = maxima - minima(i);
    matchi = find(max_diff < 0);
    if (isempty(matchi))
        max_pairs(i) = 1;
    else
        max_pairs(i) = maxima(matchi(end));
    end
end
if 0
    if (~isfield(handles, 'maxMarks'))
        handles.min_marks = line('parent', handles.im_axis, 'xdata', handles.time2(minima), 'ydata', correctedTrace(minima), ...
            'Marker', 'x', 'LineStyle', 'none', 'Color', 'r');
         handles.max_marks = line('parent', handles.im_axis, 'xdata', handles.time2(max_pairs), 'ydata', correctedTrace(max_pairs), ...
            'Marker', 'x', 'LineStyle', 'none', 'Color', 'b');
    else
        set(handles.min_marks,  'xdata', handles.time2(minima), 'ydata', correctedTrace(minima));
        set(handles.max_marks,  'xdata', handles.time2(max_pairs), 'ydata', correctedTrace(max_pairs));
    end
end

for i=1:length(minima)
    miniWindow = minima(i) + analVect;
    if (miniWindow > 0)
        if (miniWindow <= length(correctedTrace))
            miniTrace = handles.trace2_filt(miniWindow);
            miniMin = min(miniTrace(analCenter-5:analCenter+5));
            miniMax = handles.trace2_filt(max_pairs(i));
            [onsetInd, slopeTrace] = detectSingleMiniOnset(miniTrace, handles.dt2, handles.slopeThreshold);
            if ~isnan(onsetInd) %passes the slope threshold
                [onsetInd, slopeTrace] = detectSingleMiniOnset(miniTrace, handles.dt2, handles.onsetSlope);
                onsetVal = miniTrace(onsetInd);
                riseTimeSamp = handles.riseTime ./ handles.dt2;
                decayTimeSamp = handles.decayTime ./ handles.dt2;
                rise = slopeTrace(analCenter - (riseTimeSamp:-1:1) - 4);
                decay = slopeTrace(analCenter + (1:decayTimeSamp) + 4);
                if ((miniMin - miniMax) < handles.riseThreshold) %passed the rise threshold
                    if ((mean(rise < 0)>.9) && (mean(decay > 0)>.9))
                        handles.miniOnsets(handles.numMinis+1) =  analVect(onsetInd)+minima(i);
                        handles.miniPeaks(handles.numMinis+1) = minima(i);
                        handles.numMinis = handles.numMinis + 1;
                    end
                end
            end
        end
    end
    
end
handles = updateMiniTraces(hObject, handles);
guidata(hObject, handles);   

% % The dumb loop way of detecting threshold crossings, just stepping through
% % the trace point by point, with a sliding window centered on the 
% i = abs(analVect(1)) +1;
% %i=1;
% while i <= (length(correctedTrace) - windowLen/2)
%     if(correctedTrace(i) <= handles.threshold)
%         [localMin, localMinIndex] = min(correctedTrace(i+windowVect));
%         localMinIndex = windowVect(localMinIndex)+i;
%         startIndex = localMinIndex-riseSamples;
%         if startIndex < 1
%             startIndex = 1;
%         end
%         startVal = correctedTrace(startIndex);
%         difft = correctedTrace(localMinIndex) - startVal;
%         riseTimeMin = min(correctedTrace(startIndex:localMinIndex-1)); %within a rise time, it is a downward trace
%         %found an event, checking for a min deflection and that there are
%         %no intervening peaks
%         [onsetInd, slopeTrace] = detectSingleMiniOnset(handles.trace2_filt(i+analVect), handles.dt2, handles.slopeThreshold);
%         if (difft <= handles.riseThreshold && riseTimeMin > localMin && ~isnan(onsetInd)) 
%             handles.miniPeaks(handles.numMinis +1) = localMinIndex;
%             [onsetInd, slopeTrace] = detectSingleMiniOnset(handles.trace2_filt(i+analVect), handles.dt2, handles.onsetSlope);
%             handles.miniOnsets(handles.numMinis+1) =  analVect(onsetInd)+i;
%             handles.numMinis = handles.numMinis + 1;
%             %i = localMinIndex + windowLen/2; %skip to the peak
%             i = localMinIndex;
%         end
%     end
%     %i = i+windowLen/2;
%     i = i+1;
% end
% %uniquify the detected minis
% %[b, m, n] = unique(handles.miniPeaks, 'first');
% %handles.miniPeaks = b;
% %handles.numMinis = length(handles.miniPeaks);
% handles = updateMiniTraces(hObject, handles);
% guidata(hObject, handles);

% --------------------------------------------------------------
function handles = updateMiniTraces(hObject, handles)
% Function to update the miniTraces that are stored.          
windowLen = handles.miniAnalWindow/handles.dt2;
windowVect = (0:windowLen) - windowLen/2;
handles.miniSlopes = NaN*ones(length(windowVect)-1, handles.numMinis);
handles.minis = NaN*ones(length(windowVect), handles.numMinis);
handles.miniHeights = zeros(handles.numMinis,1);
for i=1:handles.numMinis
    miniWindow = windowVect+handles.miniOnsets(i);
    neg_inds = find(miniWindow < 1);
    if (~isempty(neg_inds))
        miniWindow(neg_inds) = 1;
    end
    miniTrace = handles.trace2_filt(miniWindow) - handles.trace2_filt(handles.miniOnsets(i)); %recenter the onsets
    if (~isempty(neg_inds))
        miniTrace(neg_inds) = NaN;
    end
    handles.miniSlopes(:,i) = computeMiniSlope(miniTrace, handles.dt2);
    handles.minis(:,i) = miniTrace;
    handles.miniHeights(i) = handles.trace2_filt(handles.miniPeaks(i)) - handles.trace2_filt(handles.miniOnsets(i));
end
%guidata(hObject,handles);


% --------------------------------------------------------------
function [onsetInd, slopeTrace] = detectSingleMiniOnset(miniTrace, dt, slopeThreshold)
% function [onsetInd, slopeTrace] = detectMiniOnset(miniTrace)
% This takes a mini trace and returns the onset index.  If the event
% given doesn't pass the threshold for onset slope, then the onsetInd
% returned is NaN.
edgeSize = 10;  %the number of samples to ignore on the edges
% filter, zero lag, with boxcar of filtWindSize
slopeTrace = computeMiniSlope(miniTrace, dt);
sigSlopeInds = find(slopeTrace < slopeThreshold);
if (isempty(sigSlopeInds))
    onsetInd = NaN;
    return;
else
    onsetInd = sigSlopeInds(1);
end


% -------------------------------------------------
function handles = detectMiniOnsets(handles)
%
% detection of the onsets.  Gonna try to use the derivatives to get the
% onset using a derivative threshold.
windowVect = 1:(handles.miniAnalWindow/handles.dt2);
windowVect = windowVect - floor((length(windowVect)/2));
centeri = floor((length(windowVect)/2)) + 1;
wvLen = length(windowVect);
filtWindSize = 9; %filter before taking derivative
edgeSize = 10;  %the number of samples to ignore on the edges
pb = 1;
slopeThreshold = -.015;
miniArray = NaN*zeros(handles.numMinis, wvLen*2);
figure; traceah = axes; hold on;
figure; derivah = axes; hold on;
for i = 1:handles.numMinis
    miniTrace = handles.trace2(windowVect+handles.miniPeaks(i));
    miniTraceFilt = filtfilt(ones(1,filtWindSize)/filtWindSize,1,miniTrace);
    %miniTraceFilt = miniTraceFilt(floor(filtWindSize/2):end);
    %miniTraceFilt = medfilt1(miniTrace, 3);
    mtSlope = diff(miniTraceFilt)/handles.dt2;
    sigSlopeInds = [];
    sigSlopeInds = find(mtSlope < slopeThreshold);
    if (~isempty(sigSlopeInds))
        threshCrossing = sigSlopeInds(1); %find the first point of crossing
        baseVm = mean(miniTrace(1:threshCrossing));
        miniArray(i, (1:wvLen)+wvLen-threshCrossing) = miniTrace-baseVm;
    end
    if pb
        plot(derivah, mtSlope(1+edgeSize:end-edgeSize));
        %plot(traceah, miniTrace, 'k'); %plot(traceah, miniTraceFilt, 'b');
        plot(traceah, miniArray(i,:), 'k');
        if (~isempty(sigSlopeInds))
            %plot(traceah, threshCrossing, miniTraceFilt(threshCrossing), 'rx');
        end
    end
end
meanMini = nanmean2(miniArray);
plot(traceah,meanMini, 'r', 'LineWidth', 2);


% --------------------------------------------------------------
function deleteClickedMini(src, evt, hObject)
% Deletes the clicked mini from the miniPeaks
handles = guidata(hObject);
cp = get(gca, 'CurrentPoint'); %get last clicked X position

if (~isempty(handles.miniPeaks))
    %find a close mini
    miniTimes = handles.time2_lp(handles.miniPeaks);
    timeDiff = cp(1) - miniTimes;
    [mintd minind] = min(abs(timeDiff));
    if (abs(mintd) < 10) % found the match
        newMini = [handles.miniPeaks(1:(minind-1)) handles.miniPeaks((minind+1):end)];
        handles.miniPeaks = newMini;
        newMini = [handles.miniOnsets(1:(minind-1)) handles.miniOnsets((minind+1):end)];
        handles.miniOnsets = newMini;
        handles.numMinis = handles.numMinis - 1;
    end  
end
handles = markMinis(hObject, [], handles);
guidata(hObject,handles);

    
% ----------------------------------------------------------------
function handles = markMinis(hObject, eventdata, handles)
% Creates a plot object for the current list of detected minis
%handles = guidata(hObject);
if (~isempty(handles.miniplot))
    if (~isempty(handles.miniPeaks))
        set(handles.miniplot, 'XData', handles.time2_lp(handles.miniPeaks), ...
            'YData', handles.trace2_sub(handles.miniPeaks));
    else
        set(handles.miniplot, 'xData', [], 'Ydata', []);
    end
else
    handles.miniplot =  line('Parent', handles.im_axis, ...
                             'XData', handles.time2_lp(handles.miniPeaks), ...
                             'YData', handles.trace2_sub(handles.miniPeaks), ...
                             'Color', 'b', 'Marker', 'o',  'LineStyle', 'none', ...
                             'MarkerSize', 12);
end
set(handles.miniNum_text, 'String', num2str(handles.numMinis));
set(handles.miniplot, 'ButtonDownFcn', {@deleteClickedMini, hObject}); %to be deleted on another click
%guidata(hObject,handles);


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% file = uigetfile('*.fig');
% if ~isequal(file, 0)
%     open(file);
% end
loadfiles(hObject, eventdata, handles);

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% printdlg(handles.miniAnal_fig)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% selection = questdlg(['Close ' get(handles.miniAnal_fig,'Name') '?'],...
%                      ['Close ' get(handles.miniAnal_fig,'Name') '...'],...
%                      'Yes','No','Yes');
% if strcmp(selection,'No')
%     return;
% end
% 
% delete(handles.miniAnal_fig)


function tmin_edit_Callback(hObject, eventdata, handles)
% hObject    handle to tmin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tmin_edit as text
%        str2double(get(hObject,'String')) returns contents of tmin_edit as a double
tmp = str2double(get(hObject,'String'));
if ~limitCheck(tmp, handles.xlim(2))
    disp('Max value must be greater than min value');
    return;
else
    handles.xlim(1) = tmp;
    setXLimits(hObject, eventdata, handles);
    guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function tmin_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tmin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in tmin_dec.
function tmin_dec_Callback(hObject, eventdata, handles)
% hObject    handle to tmin_dec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new = handles.xlim(1) - 50;
if limitCheck(new, handles.xlim(2))
    set(handles.tmin_edit, 'String', num2str(new));
    handles.xlim(1) = new;
    setXLimits(hObject, eventdata, handles);
    guidata(hObject, handles);
end

% --- Executes on button press in tmin_inc.
function tmin_inc_Callback(hObject, eventdata, handles)
% hObject    handle to tmin_inc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new = handles.xlim(1) + 50;
if limitCheck(new, handles.xlim(2))
    set(handles.tmin_edit, 'String', num2str(new));
    handles.xlim(1) = new;
    setXLimits(hObject, eventdata, handles);
    guidata(hObject, handles);
end



function tmax_edit_Callback(hObject, eventdata, handles)
% hObject    handle to tmax_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tmax_edit as text
%        str2double(get(hObject,'String')) returns contents of tmax_edit as
%        a double
tmp = str2double(get(hObject,'String'));
if ~limitCheck(handles.xlim(1), tmp)
    disp('Max value must be greater than min value');
    return;
else
    handles.xlim(2) = tmp;
    setXLimits(hObject, eventdata, handles);
    guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function tmax_edit_CreateFcn(hObject, eventdata, handles)


% --------------------------------------------------------------------
function display_menu_Callback(hObject, eventdata, handles)
% hObject    handle to display_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function viewOverlays_Callback(hObject, eventdata, handles)
% hObject    handle to viewOverlays (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[allMinis, allMiniSlopes] = returnAllMinis(handles);
prevNum = str2num(get(handles.prevMinis_text, 'string')); 
totalMiniNum = handles.numMinis+prevNum;
minis = zeros(size(handles.minis,1), totalMiniNum) * NaN;
windowLen = handles.miniAnalWindow/handles.dt2;
windowVect = (0:windowLen) - windowLen/2;
miniTime = ((1:length(windowVect))-1) * handles.dt2;
%handles = detectMiniOnsets(handles);

figure; hold on;
for i = 1:totalMiniNum
    neg_slopes = find(allMiniSlopes(:,i) < 0);
    %minis(i, :) = handles.trace2_sub(windowVect+handles.miniPeaks(i));
    minis(:, i) = allMinis(:, i);
    plot(miniTime, minis(:,i), 'k');
    % onset is aligned at half. Going to average till the EPSC upslope ends, then it is NaN.  
    slopeN = size(allMiniSlopes,1);
    last_half = floor(slopeN/2)+1:slopeN; 
    neg_slopes = intersect(last_half, find(allMiniSlopes(:,i) < 0));
    pos_slopes = setdiff(last_half, neg_slopes); %these are the upslopes
    jumps = find(diff(pos_slopes) > 1); %non-contiguous segments
    if isempty(jumps)
        change_ind = slopeN;
    else
        change_ind = pos_slopes(jumps(1))+1;
    end
    plot(miniTime(neg_slopes), minis(neg_slopes,i), 'b');
    plot(miniTime(pos_slopes), minis(pos_slopes,i), 'g');
    plot(miniTime(change_ind), minis(change_ind,i), 'gx');
    minis(change_ind:end, i) = NaN;
end
meanMini = nanmean2(minis, 2);
plot(miniTime, meanMini, 'r', 'Linewidth', 2);
title_str = sprintf('Mean Spontaneous Event Shape - %s', handles.expdir);
title(title_str,'Fontsize', 16);
if (handles.traceType == 1)
    ystr = '\DeltaIm (nA)';
else
    ystr = 'Inverted \DeltaVm (mV)';
end
ylabel(ystr, 'Fontsize', 14);
xlabel('Time (ms), Onset at 10ms', 'Fontsize', 14);



% -------------------- Sets the axis limits -----------
function setXLimits (hObject, eventdata, handles)

set(handles.im_axis, 'Xlim',  handles.xlim);
set(handles.vm_axis, 'Xlim', handles.xlim);
autoSetYLimits(hObject, eventdata, handles);

function autoSetYLimits(hObject, eventdata, handles)
% Just finds the min and max for each of the channels and sets the limits a
% bit beyond those
dispMin = find((handles.time1 >= handles.xlim(1)));  dispMin = dispMin(1);
dispMax = find((handles.time1 <= handles.xlim(2)));  dispMax = dispMax(end);
chanMin = min(handles.trace1(dispMin:dispMax))-.01;
chanMax = max(handles.trace1(dispMin:dispMax));
set(handles.vm_axis, 'Ylim', [chanMin chanMax]);

dispMin = find((handles.time2_lp >= handles.xlim(1)));  dispMin = dispMin(1);
dispMax = find((handles.time2_lp <= handles.xlim(2)));  dispMax = dispMax(end);
chanMin = min(handles.trace2_sub(dispMin:dispMax)) - .02;
chanMax = max(handles.trace2_sub(dispMin:dispMax)) + .02;
set(handles.im_axis, 'Ylim', [chanMin chanMax]);
%set(handles.im_axis, 'Ylim', [-4 .2]);


% ---- Checks if axis limits are valid -------
function good = limitCheck(min, max)
good = (min < max);


function decayTime_edit_Callback(hObject, eventdata, handles)
% hObject    handle to decayTime_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of decayTime_edit as text
%        str2double(get(hObject,'String')) returns contents of decayTime_edit as a double
handles.decayTime = str2double(get(hObject,'String'));
%set(handles.thresholdline, 'YData', [handles.threshold handles.threshold]);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function decayTime_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to decayTime_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function riseTime_edit_Callback(hObject, eventdata, handles)
% hObject    handle to riseTime_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of riseTime_edit as text
%        str2double(get(hObject,'String')) returns contents of riseTime_edit as a double
handles.riseTime = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function riseTime_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to riseTime_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in detect_btn.
function detect_btn_Callback(hObject, eventdata, handles)
% hObject    handle to detect_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = detectMinis2(hObject, eventdata, handles);
handles = markMinis(hObject, eventdata, handles);
guidata(hObject, handles);

% --- Executes on button press in nextTrial_btn.
function nextTrial_btn_Callback(hObject, eventdata, handles)
% hObject    handle to nextTrial_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.currn = handles.currn + 1;
handles.keepMinis = 1;
loadTrial(hObject, eventdata, handles);


% --- Executes on button press in nextTrial_btn.
function previousTrial_btn_Callback(hObject, eventdata, handles)
% hObject    handle to nextTrial_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.currn = handles.currn - 1;
if(handles.currn < 1) handles.currn = 1; end
handles.keepMinis = 1;
loadTrial(hObject, eventdata, handles);


% --- Executes on button press in nextPeriodbtn.
function nextPeriodbtn_Callback(hObject, eventdata, handles)
% hObject    handle to nextPeriodbtn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

period = handles.xlim(2) - handles.xlim(1);
handles.xlim = handles.xlim+period;
set(handles.tmin_edit, 'String', num2str(handles.xlim(1)));
set(handles.tmax_edit, 'String', num2str(handles.xlim(2)));
setXLimits(hObject, eventdata, handles);
guidata(hObject, handles);

% --- Executes on button press in previousPeriod_btn.
function previousPeriod_btn_Callback(hObject, eventdata, handles)
% hObject    handle to previousPeriod_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
period = handles.xlim(2) - handles.xlim(1);
handles.xlim = handles.xlim-period;
set(handles.tmin_edit, 'String', num2str(handles.xlim(1)));
set(handles.tmax_edit, 'String', num2str(handles.xlim(2)));
setXLimits(hObject, eventdata, handles);
guidata(hObject, handles);


% --- Executes on button press in tmax_dec.
function tmax_dec_Callback(hObject, eventdata, handles)
% hObject    handle to tmax_dec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new = handles.xlim(2) - 50;
if limitCheck(handles.xlim(1), new)
    set(handles.tmax_edit, 'String', num2str(new));
    handles.xlim(2) = new;
    setXLimits(hObject, eventdata, handles);
    guidata(hObject, handles);
end

% --- Executes on button press in tmax_inc.
function tmax_inc_Callback(hObject, eventdata, handles)
% hObject    handle to tmax_inc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new = handles.xlim(2) + 50;
if limitCheck(handles.xlim(1), new)
    set(handles.tmax_edit, 'String', num2str(new));
    handles.xlim(2) = new;
    setXLimits(hObject, eventdata, handles);
    guidata(hObject, handles);
end



function riseThreshold_edit_Callback(hObject, eventdata, handles)
% hObject    handle to riseThreshold_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of riseThreshold_edit as text
%        str2double(get(hObject,'String')) returns contents of riseThreshold_edit as a double
handles.riseThreshold = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function riseThreshold_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to riseThreshold_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function start_move(src, evt, hObject)
  set(gcbf,'WindowButtonMotionFcn','miniAnal move_thres', ...
	   'WindowButtonUpFcn', 'miniAnal stop_move');
  
% --------------------------------------------------------------------
function stop_move(src, evt, hObject)
  set(gcbf,'WindowButtonMotionFcn','', ...
	   'WindowButtonUpFcn', ''); 
  handles = guidata(gcbf); 
  thresh  = get(gco, 'Ydata');
  handles.threshold = thresh(1);
  %set(handles.decayTime_edit, 'String', num2str(handles.threshold));
  guidata(gcbf, handles);
   
  
% --------------------------------------------------------------------
function move_thres(src, evt, hObject)
  curr_pt = get(gca,'CurrentPoint');
  set(gco,'YData',[curr_pt(1,2) curr_pt(1,2)]);


% --- Executes on button press in plotVmDist_btn.
function plotVmDist_btn_Callback(hObject, eventdata, handles)
% hObject    handle to plotVmDist_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure;
ah = axes;
bin_inc = .02;
plotVectorDistribution(ah, handles.trace1, bin_inc, 'b');

% --- Executes on button press in plotImDist_btn.
function plotImDist_btn_Callback(hObject, eventdata, handles)
% hObject    handle to plotImDist_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure; 
ah = axes;
bin_inc = .01;
plotVectorDistribution(ah, handles.trace2_sub, bin_inc, 'k');


% --------------------------------------------------------------------
function viewMiniHeights_Callback(hObject, eventdata, handles)
% hObject    handle to viewMiniHeights (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[allMinis, allMiniSlopes] = returnAllMinis(handles);
prevNum = str2num(get(handles.prevMinis_text, 'string')); 
totMiniNum = prevNum + handles.numMinis;
miniHeights = zeros(totMiniNum, 1);
num = 0;
for i = 1:length(handles.prev)
    nums = num + (1:handles.prev(i).numMinis);
    miniHeights(nums) = handles.prev(i).miniHeights(1:handles.prev(i).numMinis);
    num = nums(end);
end
nums = num + (1:handles.numMinis);
miniHeights(nums) = handles.miniHeights(1:handles.numMinis);
bin_inc = .05;
figure;
ah=axes;
[bh, Y, bins] = plotVectorDistribution(ah, miniHeights, bin_inc, 'k');
medianHeight = nanmedian(miniHeights);
meanHeight = nanmean2(miniHeights);
lh = line('Parent', ah, 'Xdata', [medianHeight medianHeight], 'Ydata', [0 max(Y) + 2], 'Color', 'r', 'Linestyle', '--', 'linewidth', 2);
text(medianHeight+.1, max(Y) + 1, ['Median = ' num2str(medianHeight)]);
text(medianHeight+.1, max(Y) + 2, ['Mean = ' num2str(meanHeight)]);
figure; hold on;
plot(allMiniSlopes);


% --- Executes on selection change in traceType_menu.
function traceType_menu_Callback(hObject, eventdata, handles)
% hObject    handle to traceType_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns traceType_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from traceType_menu
tt = get(hObject, 'value');
if (tt == 2 && handles.traceType==1) %switch from current to voltage anal
    factor = handles.dv/handles.di;
elseif (tt == 1 && handles.traceType==2) %switch from voltage to current anal
    factor = handles.di/handles.dv;
else
    factor=1;
end
handles.riseThreshold = handles.riseThreshold * factor;
handles.slopeThreshold = handles.slopeThreshold * factor;
handles.threshold = handles.threshold * factor;
handles.traceType = tt; 

%take care of gui stuff
if (~isempty(handles.thresholdline)) set(handles.thresholdline, 'Ydata', [handles.threshold handles.threshold]); end
if (handles.traceType ==2) %Vm analysis
    handles.invertTrace2 = 1;
    set(handles.trace2Invert_check, 'value', get(handles.trace2Invert_check,'max'));
end
set(handles.riseThreshold_edit, 'string', num2str(handles.riseThreshold));
set(handles.slopeThreshold_edit, 'string', num2str(handles.slopeThreshold));
set(handles.decayTime_edit, 'string', num2str(handles.threshold));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function traceType_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to traceType_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

labels = {'Current trace'; 'Voltage Trace'};
set(hObject, 'string', labels);


% --- Executes on button press in trace2Invert_check.
function trace2Invert_check_Callback(hObject, eventdata, handles)
% hObject    handle to trace2Invert_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trace2Invert_check

handles.invertTrace2 = get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in trace1Invert_check.
function trace1Invert_check_Callback(hObject, eventdata, handles)
% hObject    handle to trace1Invert_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trace1Invert_check

handles.invertTrace1 = get(hObject,'Value');
guidata(hObject, handles);



function slopeThreshold_edit_Callback(hObject, eventdata, handles)
% hObject    handle to slopeThreshold_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of slopeThreshold_edit as text
%        str2double(get(hObject,'String')) returns contents of slopeThreshold_edit as a double
handles.slopeThreshold = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slopeThreshold_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slopeThreshold_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in export_btn.
function export_btn_Callback(hObject, eventdata, handles)
% hObject    handle to export_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.prev)
    prevNum = handles.currn;
else
    fn = handles.prev(1).filename;
    toks = strread(fn, '%s', 'delimiter', '/');
    fn = toks{end};
    t2 = strread(fn, '%s', 'delimiter', '_');
    prevNum = str2num(t2{2});
end
fname = sprintf('%s/%s_%i-%i_mini.mat', handles.rootdir, handles.expdir,prevNum, handles.currn);

save(fname, 'handles');



function onsetSlope_edit_Callback(hObject, eventdata, handles)
% hObject    handle to onsetSlope_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of onsetSlope_edit as text
%        str2double(get(hObject,'String')) returns contents of onsetSlope_edit as a double

handles.onsetSlope = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function onsetSlope_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to onsetSlope_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

 % --------------------------------------------------------------------
function addClickedEvent(src, evt, hObject)
  curr_pt = get(gca,'CurrentPoint');
  curr_x = curr_pt(1,1);
  handles = guidata(gcbf);
  [minval, clicked_i] = min(abs(handles.time2 - curr_x));
  windowVect = 1:(handles.miniSearchWindow/handles.dt2);
  windowVect = windowVect - floor((length(windowVect)/2));
  miniTrace = handles.trace2_filt(windowVect+clicked_i);
  [miniPeak, miniPeak_i] = min(miniTrace);
  [onsetInd, slopeTrace] = detectSingleMiniOnset(miniTrace, handles.dt2, handles.onsetSlope);
  if (isnan(onsetInd)) % could not detect onset
      %exit function without saving anything
      return;
  end
  %now let's update the structures that save the mini
  handles.numMinis = handles.numMinis+1;
  handles.miniPeaks(handles.numMinis) = windowVect(miniPeak_i) + clicked_i;
  handles.miniOnsets(handles.numMinis) = windowVect(onsetInd) + clicked_i;

  handles = markMinis(gcbf, [], handles);
  handles = updateMiniTraces(gcbf, handles);
  guidata(gcbf, handles);
 
  


% --- Executes on selection change in format_menu.
function format_menu_Callback(hObject, eventdata, handles)
% hObject    handle to format_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns format_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from format_menu
contents = get(hObject,'String')
format = contents{get(hObject, 'value')};
handles.format = format;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function format_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to format_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

labels = {'int16'; 'double'};
value = 1;
set(hObject, 'string', labels, 'value', value);





    



function trace2amp_edit_Callback(hObject, eventdata, handles)
% hObject    handle to trace2amp_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trace2amp_edit as text
%        str2double(get(hObject,'String')) returns contents of trace2amp_edit as a double
handles.trace2_amp = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function trace2amp_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trace2amp_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


