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

% Last Modified by GUIDE v2.5 01-Jul-2013 14:50:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @movieAnal_OpeningFcn, ...
                   'gui_OutputFcn',  @movieAnal_OutputFcn, ...
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

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes movieAnal wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = movieAnal_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


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

%show the average image from the movie
axes(handles.im_ax);
imshow(handles.mt.avgFrame, 'InitialMagnification', 'fit');


% Update handles structure
guidata(hObject, handles);
