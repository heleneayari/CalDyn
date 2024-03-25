function varargout = Parameters(varargin)
% PARAMETERS MATLAB code for Parameters.fig
%      PARAMETERS, by itself, creates a new PARAMETERS or raises the existing
%      singleton*.
%
%      H = PARAMETERS returns the handle to a new PARAMETERS or the handle to
%      the existing singleton*.
%
%      PARAMETERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARAMETERS.M with the given input arguments.
%
%      PARAMETERS('Property','Value',...) creates a new PARAMETERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Parameters_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Parameters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Parameters

% Last Modified by GUIDE v2.5 25-Mar-2024 15:51:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Parameters_OpeningFcn, ...
                   'gui_OutputFcn',  @Parameters_OutputFcn, ...
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


% --- Executes just before Parameters is made visible.
function Parameters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Parameters (see VARARGIN)

% Choose default command line output for Parameters

handles.output = hObject;
p = inputParser;
addRequired(p,'list_allparam');
addRequired(p, 'list_param_num'); 
addOptional(p, 'list_calc',[1 1 1]); 

parse(p, varargin{:});
handles.list_allparam=p.Results.list_allparam;
handles.list_param_num=p.Results.list_param_num;
handles.list_param_name=handles.list_allparam(handles.list_param_num(:));

handles.list_calc=p.Results.list_calc;
% get(handles.N_pks,'Value')
% set(handles.N_pks,'Value',0)

if handles.list_calc(1)
   set(handles.calc_mean,'Value',1) ;
else
    set(handles.calc_mean,'Value',0) ;
end
if handles.list_calc(2)
   set(handles.calc_median,'Value',1) ;
else
     set(handles.calc_median,'Value',0) ;
end
if handles.list_calc(3)
   set(handles.calc_std,'Value',1) ;
else
     set(handles.calc_std,'Value',0) ;
end

for ii=1:length(handles.list_param_name)
   
   eval(['set(handles.',handles.list_param_name{ii},',',char(39),'Value',char(39),',1)']);
   handles.list_param_num(strcmp(handles.list_allparam,handles.list_param_name(ii)))=1;
    
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Parameters wait for user response (see UIRESUME)
 uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Parameters_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.list_param_num;
varargout{2} = handles.list_calc;
if (isfield(handles,'closeFigure') && handles.closeFigure)
    figure1_CloseRequestFcn(hObject, eventdata, handles)
end
function figure1_CloseRequestFcn(hObject, ~, ~)
delete(hObject);



% --- Executes on button press in N_pks.
function N_pks_Callback(hObject, eventdata, handles)
% hObject    handle to N_pks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);

% --- Executes on button press in Period.
function Period_Callback(hObject, eventdata, handles)
% hObject    handle to Period (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Period
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);

% --- Executes on button press in Asc_time.
function Asc_time_Callback(hObject, eventdata, handles)
% hObject    handle to Asc_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Asc_time
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);

% --- Executes on button press in Decay_time.
function Decay_time_Callback(hObject, eventdata, handles)
% hObject    handle to Decay_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Decay_time
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);

% --- Executes on button press in Decay_time_95.
function Decay_time_95_Callback(hObject, eventdata, handles)
% hObject    handle to Decay_time_95 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Decay_time_95
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);

% --- Executes on button press in AUC.
function AUC_Callback(hObject, eventdata, handles)
% hObject    handle to AUC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AUC
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);

% --- Executes on button press in Asc_slope.
function Asc_slope_Callback(hObject, eventdata, handles)
% hObject    handle to Asc_slope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Asc_slope
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);

% --- Executes on button press in Decay_slope_0_50.
function Decay_slope_0_50_Callback(hObject, eventdata, handles)
% hObject    handle to Decay_slope_0_50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Decay_slope_0_50
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);

% --- Executes on button press in Amp_asc.
function Amp_asc_Callback(hObject, eventdata, handles)
% hObject    handle to Amp_asc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Amp_asc
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);

% --- Executes on button press in calc_mean.
function calc_mean_Callback(hObject, eventdata, handles)
% hObject    handle to calc_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of calc_mean
val=get(hObject,'Value');
handles.list_calc(1)=val;
guidata(hObject, handles);

% --- Executes on button press in Taud.
function Taud_Callback(hObject, eventdata, handles)
% hObject    handle to Taud (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Taud
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);

% --- Executes on button press in f_smpks.
function f_smpks_Callback(hObject, eventdata, handles)
% hObject    handle to f_smpks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of f_smpks
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);

% --- Executes on button press in f_multipks.
function f_multipks_Callback(hObject, eventdata, handles)
% hObject    handle to f_multipks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of f_multipks
val=get(hObject,'Value');
pos=strcmp(handles.list_allparam,'N_pks');
handles.list_param_num(pos)=val;
guidata(hObject, handles);

% --- Executes on button press in FP_duration.
function FP_duration_Callback(hObject, eventdata, handles)
% hObject    handle to FP_duration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FP_duration
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);



% --- Executes on button press in Decay_slope.
function Decay_slope_Callback(hObject, eventdata, handles)
% hObject    handle to Decay_slope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Decay_slope
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);


% --- Executes on button press in Decay_slope_50_100.
function Decay_slope_50_100_Callback(hObject, eventdata, handles)
% hObject    handle to Decay_slope_50_100 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Decay_slope_50_100
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);


% --- Executes on button press in Maxima.
function Maxima_Callback(hObject, eventdata, handles)
% hObject    handle to Maxima (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Maxima
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);


% --- Executes on button press in calc_median.
function calc_median_Callback(hObject, eventdata, handles)
% hObject    handle to calc_median (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of calc_median
val=get(hObject,'Value');
handles.list_calc(2)=val;
guidata(hObject, handles);

% --- Executes on button press in Baz_taud.
function Baz_taud_Callback(hObject, eventdata, handles)
% hObject    handle to Baz_taud (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Baz_taud
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);


% --- Executes on button press in f_medpks.
function f_medpks_Callback(hObject, eventdata, handles)
% hObject    handle to f_medpks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of f_medpks
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);




% --- Executes on button press in FP_Amp.
function FP_Amp_Callback(hObject, eventdata, handles)
% hObject    handle to FP_Amp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of FP_Amp
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);






% --- Executes on button press in calc_std.
function calc_std_Callback(hObject, eventdata, handles)
% hObject    handle to calc_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of calc_std
val=get(hObject,'Value');
handles.list_calc(3)=val;
guidata(hObject, handles);


% --- Executes on button press in Minima.
function Minima_Callback(hObject, eventdata, handles)
% hObject    handle to Minima (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Minima
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);


% --- Executes on button press in Amp_decay.
function Amp_decay_Callback(hObject, eventdata, handles)
% hObject    handle to Amp_decay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Amp_decay
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);




% --- Executes on button press in Decay_time_70.
function Decay_time_70_Callback(hObject, eventdata, handles)
% hObject    handle to Decay_time_70 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Decay_time_70
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);


% --- Executes on button press in Decay_time_50.
function Decay_time_50_Callback(hObject, eventdata, handles)
% hObject    handle to Decay_time_50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Decay_time_50
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);


% --- Executes on button press in Decay_time_30.
function Decay_time_30_Callback(hObject, eventdata, handles)
% hObject    handle to Decay_time_30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Decay_time_30
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);


% --- Executes on button press in Decay_time_20.
function Decay_time_20_Callback(hObject, eventdata, handles)
% hObject    handle to Decay_time_20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Decay_time_20
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);

% --- Executes on button press in Sig_noise.
function Sig_noise_Callback(hObject, eventdata, handles)
% hObject    handle to Sig_noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Sig_noise
val=get(hObject,'Value');
pos=strcmp(handles.list_allparam,'Sig_noise');
handles.list_param_num(pos)=val;
guidata(hObject, handles);

% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handl
handles.closeFigure = true;
uiresume(handles.figure1);
guidata(hObject, handles);


% --- Executes on button press in Std_decay_slope.
function Std_decay_slope_Callback(hObject, eventdata, handles)
% hObject    handle to Std_decay_slope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Std_decay_slope
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);


% --- Executes on button press in Decay_time_90.
function Decay_time_90_Callback(hObject, eventdata, handles)
% hObject    handle to Decay_time_90 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Decay_time_90
val=get(hObject,'Value');
nom=get(hObject,'Tag');
pos=strcmp(handles.list_allparam,nom);
handles.list_param_num(pos)=val;
guidata(hObject, handles);
