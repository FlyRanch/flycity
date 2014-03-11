function varargout = Fly_cockpit(varargin)
% FLY_COCKPIT MATLAB code for Fly_cockpit.fig
%      FLY_COCKPIT, by itself, creates a new FLY_COCKPIT or raises the existing
%      singleton*.
%
%      H = FLY_COCKPIT returns the handle to a new FLY_COCKPIT or the handle to
%      the existing singleton*.
%
%      FLY_COCKPIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FLY_COCKPIT.M with the given input arguments.
%
%      FLY_COCKPIT('Property','Value',...) creates a new FLY_COCKPIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Fly_cockpit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Fly_cockpit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Fly_cockpit

% Last Modified by GUIDE v2.5 24-Jan-2013 15:37:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Fly_cockpit_OpeningFcn, ...
                   'gui_OutputFcn',  @Fly_cockpit_OutputFcn, ...
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


% --- Executes just before Fly_cockpit is made visible.
function Fly_cockpit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Fly_cockpit (see VARARGIN)

% Choose default command line output for Fly_cockpit
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Fly_cockpit wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Fly_cockpit_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in startbutton.
function startbutton_Callback(hObject, eventdata, handles)
% hObject    handle to startbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
