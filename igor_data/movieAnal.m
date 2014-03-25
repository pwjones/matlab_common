function varargout = movieAnal(varargin)
% MOVIEANAL MATLAB code for movieAnal.fig
%      MOVIEANAL, by itself, creates a new MOVIEANAL or raises the existing
%      singleton*.
%
%      H = MOVIEANAL returns the handle to a new MOVIEANAL or the handle to
%      the existing singleton*.
%
%      MOVIEANAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOVIEANAL.M with the given input arguments.
%
%      MOVIEANAL('Property','Value',...) creates a new MOVIEANAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before movieAnal_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to movieAnal_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help movieAnal

% Last Modified by GUIDE v2.5 03-Jul-2013 14:33:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @movieAnal_OpeningFcn, ...
                   'gui_OutputFcn',  @movieAnal_OutputFcn, ...
                   'gui_CloseRequestFcn', @movieAnal_CloseRequestFcn, ...
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


% --- Executes just before movieAnal is made visible.
function movieAnal_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to movieAnal (see VARARGIN)

% Choose default command line output for movieAnal
handles.output = hObject;

handles.mt = []; % the MouseTracker object
handles.frameRange = []; %the frame range that we want to work with
handles.selRect = []; % the coordinates [x y w h] of the selection rectangle
handles.selRecth = []; % handle to area rectangle object
handles.areaFrames = []; %the frames where the position is within selected area
handles.selFrames = []; % the intersection of areaFrames and frameRange
handles.framesInMovie = [];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes movieAnal wait for user response (see UIRESUME)
%setappdata(hObject, 'waiting', 1);
%uiwait(hObject);

% Update handles structure
%guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = movieAnal_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.selFrames;
varargout{2} = 'Why doesnt this save?';

function movieAnal_CloseRequestFcn(hObject,eventdata,handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%delete(hObject);
%Check appdata flag to see if the main GUI is in a wait state

% --- Executes on button press in load_btn.
function load_btn_Callback(hObject, eventdata, handles)
% hObject    handle to load_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load the experiment file that we want
movie_folder = '~pwjones/Movies/mouse_training/';
[filename, movie_folder] = uigetfile([movie_folder '*.*']);
if filename == 0 % the user has canceled the file selection
    return;
end
handles.mt = MouseTrackerUnder2([movie_folder, filename]);
handles.frameRange = 1:handles.mt.nFrames;
handles.framesInMovie = handles.frameRange(end);
set(handles.firstFrame_edit, 'String', num2str(handles.frameRange(1)));
set(handles.lastFrame_edit, 'String', num2str(handles.frameRange(end)));


%show the average image from the movie
axes(handles.im_ax);
imshow(handles.mt.avgFrame, 'InitialMagnification', 'fit');

handles.selRect = [1 1 handles.mt.width handles.mt.height];
handles.selFrames = handles.frameRange;
handles.areaFrames = handles.selFrames;


% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in plotPos_btn.
function plotPos_btn_Callback(hObject, eventdata, handles)
% hObject    handle to plotPos_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.im_ax);
if ~isempty(handles.mt)
    handles.selFrames = intersect(handles.areaFrames, handles.frameRange);
    if ~isempty(handles.selFrames)
        handles.mt.plotPosition(handles.selFrames, handles.im_ax, 1);
        handles.selRecth = rectangle('Position', handles.selRect, 'LineWidth',2 ,'LineStyle','--', 'EdgeColor','w');   
    else
        handles.mt.plotPosition(1, handles.im_ax, 1);%just give one frame, at most a single dot.
        handles.selRecth = rectangle('Position', handles.selRect, 'LineWidth',2 ,'LineStyle','--', 'EdgeColor','w');
    end
end
% Update handles structure
guidata(hObject, handles);


function firstFrame_edit_Callback(hObject, eventdata, handles)
% hObject    handle to firstFrame_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of firstFrame_edit as text
%        str2double(get(hObject,'String')) returns contents of firstFrame_edit as a double
f = checkFrameLimits(handles, str2double(get(hObject,'String')));
fe = handles.frameRange(end);
handles.frameRange = [f:fe];

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function firstFrame_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to firstFrame_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lastFrame_edit_Callback(hObject, eventdata, handles)
% hObject    handle to lastFrame_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lastFrame_edit as text
%        str2double(get(hObject,'String')) returns contents of lastFrame_edit as a double
f = checkFrameLimits(handles, str2double(get(hObject,'String')));
fb = handles.frameRange(1);
handles.frameRange = [fb:f];


% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function lastFrame_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lastFrame_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function im_ax_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to im_ax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get initial button press point
% p1 = get(hObject, 'CurrentPoint');
% rect = [p1(1,1) p1(1,2) 50 100]
% [r2] = dragrect(rect);
% disp(num2str(r2));


% --- Executes on button press in selectArea_btn.
function selectArea_btn_Callback(hObject, eventdata, handles)
% hObject    handle to selectArea_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set(handles.im_ax, 'ButtonDownFcn', 'im_ax_ButtonDownFcn');
% k=waitforbuttonpress;
% if ~k %mouse button pressed
%     p1 = get(hObject, 'CurrentPoint');
%     rect = [p1(1,1) p1(1,2) 50 100]
%     [r2] = dragrect(rect);  
%     disp(num2str(r2));  
% end
axes(handles.im_ax);
view(2);
w = waitforbuttonpress;
if w == 0 % mouse click
    p1 = get(handles.im_ax, 'CurrentPoint');
    finalRect = rbbox;
    p2 = get(handles.im_ax,'CurrentPoint'); 
    p1 = [p1(1,1) p1(1,2)];
    p2 = [p2(1,1) p2(1,2)];
    if p1 == p2 %allow the user to click on the second point if the drag didn't register-rbbox sometimes sucks
        disp('Area selected has no pixels');
        w = waitforbuttonpress;
        p2 = get(handles.im_ax,'CurrentPoint'); 
        p2 = [p2(1,1) p2(1,2)];
    end
    if p1 ~= p2    
        orig = min(p1,p2);             % calculate locations
        offset = abs(p1-p2);
        rect = [orig, offset];
        handles.selRect = rect;
        if(~isempty(handles.selRecth) && ishandle(handles.selRecth))
            delete(handles.selRecth);
        end
        handles.selRecth = rectangle('Position', rect, 'LineWidth',2 ,'LineStyle','--', 'EdgeColor','w');    
    end
else % if keyboard input, clear the area
    if ~isempty(handles.mt)
        handles.selRect = [1 1 handles.mt.width handles.mt.height];
    end
    if(~isempty(handles.selRecth) && ishandle(handles.selRecth))
        delete(handles.selRecth);
    end
end
r = handles.selRect;
fi = find(handles.mt.nosePos(:,1) >= r(1) & handles.mt.nosePos(:,2) >= r(2) & handles.mt.nosePos(:,1) <= r(1)+r(3) & handles.mt.nosePos(:,2) <= r(2)+r(4));
handles.areaFrames = fi;
fi = intersect(fi, handles.frameRange);
handles.selFrames = fi;

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in returnFrames_btn.
function varargout = returnFrames_btn_Callback(hObject, eventdata, handles)
% hObject    handle to returnFrames_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
varargout{1} = handles.selFrames;


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if getappdata(handles.figure1,'waiting')
    % The GUI is still in UIWAIT, so call UIRESUME and return
    uiresume(hObject);
    setappdata(handles.figure1,'waiting',0)
else
    % The GUI is no longer waiting, so destroy it now.
    delete(hObject);
end

function retNum = checkFrameLimits(handles, num)
% Utility function to fix out of limit settings
retNum = num;
if retNum < 1
    retNum = 1;
elseif retNum > handles.framesInMovie
    retNum = handles.framesInMovie;
end
    
