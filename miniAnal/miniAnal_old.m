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

% Last Modified by GUIDE v2.5 27-Oct-2006 18:26:16

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
handles.di = 10/(2^16)*.5; % (Volts/bit) * .5nA/Volt sensitivity;
handles.filterWidth = 10; %number of samples to include in the median filter
%filter kernel for the low pass filtering to get a running baseline 
handles.stdg = 400; %400 samples is a gaussian with std 20ms
filtx = -3*handles.stdg:1:3*handles.stdg;
filty = my_normpdf(filtx,0,handles.stdg);
handles.filty = filty/sum(filty);
handles.riseTime = 5; %ms, event must rise in less time than this
handles.threshold = -.05; %nA, event must be more negative than this
handles.miniWindow = 20;
handles.miniList = [];
handles.numMinis = 0;
% handles to different plot items
handles.vmplot = [];
handles.implot = [];
handles.implot_filt = [];
 handles.miniplot = [];
 handles.thresholdplot = [];
handles.rootdir = [];
handles.expdir = [];
handles.fname = [];
handles.currn = [];
handles.xlim = [0 500];
% need to divide really long traces into segments and only plot a section
% in order to speed up the gui
handles.plotSection = 1;
handles.MAX_PLOT_LENGTH = 1e5;
handles.MAX_TRACE_TIME = 3500;

%set up the limit windows for the display
set(handles.tmin_edit, 'String', num2str(handles.xlim(1)));
set(handles.tmax_edit, 'String', num2str(handles.xlim(2)));
set(handles.riseTime_edit, 'String', num2str(handles.riseTime));
set(handles.threshold_edit, 'String', num2str(handles.threshold));

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
  def_rootdir = '~pwjones/glab/data/lgmd/';
end;

if ( ~isempty(handles.expdir) )
  def_expdir = handles.expdir;
else 
  def_expdir = '';
end;

if ( ~isempty(handles.fname) )
  def_fname = handles.fname;
else
  def_fname = 'led';
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

loadTrial(hObject, eventdata, handles);


% -------------------------------------------------------------
 function loadTrial(hObject, eventdata, handles)
% Loads the trial specified in the handles.structure

% clear any previous data
handles.trace1 = [];
handles.trace2 = [];
handles.trace2_filt = [];
handles.trace2_lp = [];
handles.time1 = [];
handles.time2 = [];
handles.time2_lp = [];

%load the channel 1 file
chan1file = sprintf('%s/%s/%s_%i_chan_2',handles.rootdir,handles.expdir, ...
		    handles.fname, handles.currn);
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
handles.trace1 = fread(fid,inf,'int16');
fclose(fid);
handles.trace1 = handles.trace1 * handles.dv; %convert A/D values to mV

%load the channel 2 file
chan2file = sprintf('%s/%s/%s_%i_chan_3',handles.rootdir,handles.expdir, ...
        handles.fname, handles.currn);
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
handles.trace2 = fread(fid, inf, 'int16');
fclose(fid);
handles.trace2 = handles.trace2 * handles.di; %convert A/D values to nA
handles.trace2 = handles.trace2;
handles.trace2_filt = medfilt2(handles.trace2, [handles.filterWidth, 1]); %filter the trace2
%come up with the proper time vectors
handles.time1 = (0:(length(handles.trace1)-1)) * handles.dt1;
handles.time2 = (0:(length(handles.trace2)-1)) * handles.dt2;
% The running avg time vector must be shorter by the half width of the
% filter
handles.time2_lp = handles.time2(1:(length(handles.time2)-3*handles.stdg));
%gaussian smoothing
filtbuff = fftfilt2(handles.filty,handles.trace2);
%correct for the time shift generated by the smoothing operation
handles.trace2_lp = filtbuff((1+3*handles.stdg):end);
% also generate a running baseline subtracted trace
handles.trace2_sub = handles.trace2_filt(1:length(handles.time2_lp)) - handles.trace2_lp;

%come up with the proper time vectors
handles.time1 = (0:(length(handles.trace1)-1)) * handles.dt1;
handles.time2 = (0:(length(handles.trace2)-1)) * handles.dt2;

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
                             'YData', handles.trace1, ...
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
    set(handles.thresholdplot, 'XData', [0 handles.MAX_TRACE_TIME], ...
        'YData', [handles.threshold handles.threshold]);
else
%     handles.implot = line('Parent', handles.im_axis, ...
%                              'XData', handles.time2, ...
%                              'YData', handles.trace2, ...
%                              'Color', 'k');
    handles.implot_filt = line('Parent', handles.im_axis, ...
                             'XData', handles.time2_lp, ...
                             'YData', handles.trace2_sub, ...
                             'Color', 'k');
    handles.implot_lp = line('Parent', handles.im_axis, ...
                             'XData', handles.time2_lp, ...
                             'YData', handles.trace2_lp, ...
                             'Color', [.5 .5 .5]);
     handles.thresholdplot = line('Parent', handles.im_axis, ...
                              'XData', [0 handles.MAX_TRACE_TIME], ...
                             'YData', [handles.threshold handles.threshold], ...
                             'Color', 'r');
end; 
set(handles.currn_text, 'String', num2str(handles.currn));
set(handles.im_axis, 'Xlim',handles.xlim, 'Ylim', [-Inf Inf]);
autoSetYLimits(hObject, eventdata, handles);

guidata(hObject,handles); %saves the altered handles data structure for fig h.

% --------------------------------------------------------------------
function detectMinis(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Goes through the trace2 and detects the miniature events 

handles.miniList = [];
handles.numMinis = 0;
correctedTrace = handles.trace2_sub;
riseSamples = handles.riseTime / handles.dt2;
windowLen = handles.miniWindow/handles.dt2;
windowVect = (0:windowLen) - windowLen/2;
% The dumb loop way of detecting threshold crossings
i = abs(windowVect(1)) +1;
while i <= (length(correctedTrace) - windowLen/2)
    if(correctedTrace(i) <= handles.threshold)
        [localMin, localMinIndex] = min(correctedTrace(i+windowVect));
        localMinIndex = windowVect(localMinIndex)+i;
        startIndex = localMinIndex-riseSamples;
        if startIndex < 1
            startIndex = 1;
        end
        startVal = correctedTrace(startIndex);
        diff = correctedTrace(localMinIndex) - startVal;
        if diff <= handles.threshold %found an event
            handles.miniList(handles.numMinis +1) = localMinIndex;
            handles.numMinis = handles.numMinis + 1;
        end
        if i < localMinIndex
            i = localMinIndex; %skip to the peak
        end
    end
    i = i+1;
end
%uniquify the detected minis
[b, m, n] = unique(handles.miniList, 'first');
handles.miniList = b;
handles.numMinis = length(handles.miniList);
guidata(hObject,handles);

% --------------------------------------------------------------
function deleteClickedMini(src, evt, hObject)
% Deletes the clicked mini from the miniList
handles = guidata(hObject);
cp = get(gca, 'CurrentPoint');

if (~isempty(handles.miniList))
    %find a close mini
    for i=1:length(handles.miniList)
        miniTime = handles.time2_lp(handles.miniList(i));
        timeDiff = cp(1) - miniTime;
        if (abs(timeDiff) < 10) % found the match
            newMini = [handles.miniList(1:(i-1)) handles.miniList((i+1):end)];
            handles.miniList = newMini;
            handles.numMinis = handles.numMinis - 1;
            break;
        end
    end
end
guidata(hObject,handles);
markMinis(hObject, [], handles);
    
% ----------------------------------------------------------------
function markMinis(hObject, eventdata, handles)
% Creates a plot object for the current list of detected minis
handles = guidata(hObject);
if ( ~isempty(handles.miniplot))
    set(handles.miniplot, 'XData', handles.time2_lp(handles.miniList), ...
            'YData', handles.trace2_sub(handles.miniList));
else
    handles.miniplot =  line('Parent', handles.im_axis, ...
                             'XData', handles.time2_lp(handles.miniList), ...
                             'YData', handles.trace2_sub(handles.miniList), ...
                             'Color', 'b', 'Marker', 'o',  'LineStyle', 'none', ...
                             'MarkerSize', 12);
end
set(handles.miniNum_text, 'String', num2str(handles.numMinis));
set(handles.miniplot, 'ButtonDownFcn', {@deleteClickedMini, hObject}); %to be deleted on another click
guidata(hObject,handles);

%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------- Sets the axis limits ----------------
%%%%%%%%%%%%%%%%%%%%%%%%
function setXLimits (hObject, eventdata, handles)

set(handles.im_axis, 'Xlim',  handles.xlim);
set(handles.vm_axis, 'Xlim', handles.xlim);
autoSetYLimits(hObject, eventdata, handles)

% --------------------------------------------------------------------
function autoSetYLimits(hObject, eventdata, handles)
% Just finds the min and max for each of the channels and sets the limits a
% bit beyond those
dispMin = find((handles.time1 >= handles.xlim(1)));  dispMin = dispMin(1);
dispMax = find((handles.time1 <= handles.xlim(2)));  dispMax = dispMax(end);
chanMin = min(handles.trace1(dispMin:dispMax));
chanMax = max(handles.trace1(dispMin:dispMax));
set(handles.vm_axis, 'Ylim', [chanMin chanMax]);

dispMin = find((handles.time2_lp >= handles.xlim(1)));  dispMin = dispMin(1);
dispMax = find((handles.time2_lp <= handles.xlim(2)));  dispMax = dispMax(end);
chanMin = min(handles.trace2_sub(dispMin:dispMax)) - .02;
chanMax = max(handles.trace2_sub(dispMin:dispMax)) + .02;
set(handles.im_axis, 'Ylim', [chanMin chanMax]);

% ---- Checks if axis limits are valid -------
function good = limitCheck(min, max)
good = (min < max);



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
% hObject    handle to tmax_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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


function threshold_edit_Callback(hObject, eventdata, handles)
% hObject    handle to threshold_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshold_edit as text
%        str2double(get(hObject,'String')) returns contents of threshold_edit as a double
handles.threshold = str2double(get(hObject,'String'));
set(handles.thresholdplot, 'YData', [handles.threshold handles.threshold]);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function threshold_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold_edit (see GCBO)
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
detectMinis(hObject, eventdata, handles);
markMinis(hObject, eventdata, handles);


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

minis = zeros(handles.numMinis, handles.miniWindow/handles.dt2) * NaN;
windowVect = 1:(handles.miniWindow/handles.dt2);
windowVect = windowVect - floor((length(windowVect)/2));
miniTime = ((1:length(windowVect))-1) * handles.dt2;

figure; hold on;

for i = 1:handles.numMinis
    minis(i, :) = handles.trace2_sub(windowVect+handles.miniList(i));
    plot(miniTime, minis(i,:), 'k');
end
meanMini = nanmean2(minis);
plot(miniTime, meanMini, 'r');


% --- Executes on button press in nextTrial_btn.
function nextTrial_btn_Callback(hObject, eventdata, handles)
% hObject    handle to nextTrial_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.currn = handles.currn + 1;
loadTrial(hObject, eventdata, handles);



