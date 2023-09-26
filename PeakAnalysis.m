function varargout = PeakAnalysis(varargin)
% PEAKANALYSIS MATLAB code for PeakAnalysis.fig
%      PEAKANALYSIS, by itself, creates a new PEAKANALYSIS or raises the existing
%      singleton*.
%
%      H = PEAKANALYSIS returns the handle to a new PEAKANALYSIS or the handle to
%      the existing singleton*.
%
%      PEAKANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PEAKANALYSIS.M with the given input arguments.
%
%      PEAKANALYSIS('Property','Value',...) creates a new PEAKANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PeakAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PeakAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PeakAnalysis

% Last Modified by GUIDE v2.5 24-Sep-2023 18:34:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PeakAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @PeakAnalysis_OutputFcn, ...
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

function PeakAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
set(handles.figure1,'toolbar','figure');

p = inputParser;
addRequired(p, 'PK'); %should be required..
addRequired(p, 'i'); %should be required..
addRequired(p, 'ResDir'); %should be required..
parse(p, varargin{:});
handles.PK = p.Results.PK;
handles.i = p.Results.i;
handles.ResDir=p.Results.ResDir;
handles.col=[0.5 0.5 0.5];

set(handles.pol_order,'string',num2str(handles.PK.vector_filtering_polynomial_order(handles.i)))
set(handles.frame_length,'string',num2str(handles.PK.vector_filtering_frame_length(handles.i)))
set(handles.smooth_length,'string',num2str(handles.PK.sm(handles.i)))
set(handles.prop,'string',num2str(handles.PK.prop(handles.i)))



handles.PK.Filter(handles.i);
handles.PK.CalculateParameters(handles.i);
% handles.PK.Thresh(handles.i);
% handles.PK.Calculate(handles.i);

pos=handles.PK.posper(:,:,handles.i);
M=handles.PK.Mper(:,:,handles.i);
hold(handles.axes_image,'on');
plot(handles.PK.matrix_rough_fluorescences(:,handles.i),'color',handles.col,'parent',handles.axes_image)
plot(handles.PK.matrix_filtered_fluorescences(:,handles.i),'linewidth',2,'color','b','parent',handles.axes_image)
plot(handles.PK.xmc(:,handles.i),handles.PK.mmvg(:,handles.i),'+c','parent',handles.axes_image)
plot(handles.PK.xMc(:,handles.i),handles.PK.M(:,handles.i),'+y','parent',handles.axes_image)
plot(handles.PK.xmr(:,handles.i),handles.PK.mmvd(:,handles.i),'+g','parent',handles.axes_image)
plot(handles.PK.xMr(:,handles.i),handles.PK.M(:,handles.i),'+r','parent',handles.axes_image)
plot(handles.PK.posm(:,handles.i),handles.PK.ms(:,handles.i),'+k','parent',handles.axes_image)
plot(handles.PK.posM(:,handles.i),handles.PK.Ms(:,handles.i),'+k','parent',handles.axes_image)
plot(pos(:),M(:),'+m','parent',handles.axes_image)


guidata(hObject, handles);
uiwait(handles.figure1);

function varargout = PeakAnalysis_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.PK;
if (isfield(handles,'closeFigure') && handles.closeFigure)
    figure1_CloseRequestFcn(hObject, eventdata, handles)
end
function figure1_CloseRequestFcn(hObject, ~, ~)
delete(hObject);
function Discard_Callback(hObject, eventdata, handles)

guidata(hObject, handles);

function pushbutton_save_Callback(hObject, eventdata, handles)
hgexport(handles.figure1, [handles.ResDir, 'Peak_Analysis_cell', num2str(handles.i),'.png'], hgexport('factorystyle'), 'Format', 'png');


handles.closeFigure = true;
uiresume(handles.figure1);
guidata(hObject, handles);





function pol_order_Callback(hObject, eventdata, handles)
% hObject    handle to pol_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pol_order as text
%        str2double(get(hObject,'String')) returns contents of pol_order as a double
input = str2double(get(hObject,'string'));
handles.PK.vector_filtering_polynomial_order(handles.i:end)=input;
handles.PK.Filter(handles.i);
handles.PK.CalculateParameters(handles.i);
% handles.PK.Thresh(handles.i);
% handles.PK.Calculate(handles.i);
cla(handles.axes_image)
pos=handles.PK.posper(:,:,handles.i);
M=handles.PK.Mper(:,:,handles.i);
hold(handles.axes_image,'on');
plot(handles.PK.matrix_rough_fluorescences(:,handles.i),'color',handles.col,'parent',handles.axes_image)
plot(handles.PK.matrix_filtered_fluorescences(:,handles.i),'linewidth',2,'color','b','parent',handles.axes_image)
plot(handles.PK.xmc(:,handles.i),handles.PK.mmvg(:,handles.i),'+c','parent',handles.axes_image)
plot(handles.PK.xMc(:,handles.i),handles.PK.M(:,handles.i),'+y','parent',handles.axes_image)
plot(handles.PK.xmr(:,handles.i),handles.PK.mmvd(:,handles.i),'+g','parent',handles.axes_image)
plot(handles.PK.xMr(:,handles.i),handles.PK.M(:,handles.i),'+r','parent',handles.axes_image)
plot(handles.PK.posm(:,handles.i),handles.PK.ms(:,handles.i),'+k','parent',handles.axes_image)
plot(handles.PK.posM(:,handles.i),handles.PK.Ms(:,handles.i),'+k','parent',handles.axes_image)
plot(pos(:),M(:),'+m','parent',handles.axes_image)
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pol_order_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pol_order (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frame_length_Callback(hObject, eventdata, handles)
% hObject    handle to frame_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_length as text
%        str2double(get(hObject,'String')) returns contents of frame_length as a double
input = str2double(get(hObject,'string'));
handles.PK.vector_filtering_frame_length(handles.i:end)=input;
handles.PK.Filter(handles.i);
handles.PK.CalculateParameters(handles.i);
% handles.PK.Thresh(handles.i);
% handles.PK.Calculate(handles.i);
cla(handles.axes_image)
pos=handles.PK.posper(:,:,handles.i);
M=handles.PK.Mper(:,:,handles.i);
hold(handles.axes_image,'on');
plot(handles.PK.matrix_rough_fluorescences(:,handles.i),'color',handles.col,'parent',handles.axes_image)
plot(handles.PK.matrix_filtered_fluorescences(:,handles.i),'linewidth',2,'color','b','parent',handles.axes_image)
plot(handles.PK.xmc(:,handles.i),handles.PK.mmvg(:,handles.i),'+c','parent',handles.axes_image)
plot(handles.PK.xMc(:,handles.i),handles.PK.M(:,handles.i),'+y','parent',handles.axes_image)
plot(handles.PK.xmr(:,handles.i),handles.PK.mmvd(:,handles.i),'+g','parent',handles.axes_image)
plot(handles.PK.xMr(:,handles.i),handles.PK.M(:,handles.i),'+r','parent',handles.axes_image)
plot(handles.PK.posm(:,handles.i),handles.PK.ms(:,handles.i),'+k','parent',handles.axes_image)
plot(handles.PK.posM(:,handles.i),handles.PK.Ms(:,handles.i),'+k','parent',handles.axes_image)
plot(pos(:),M(:),'+m','parent',handles.axes_image)
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function frame_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function smooth_length_Callback(hObject, eventdata, handles)
% hObject    handle to smooth_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of smooth_length as text
%        str2double(get(hObject,'String')) returns contents of smooth_length as a double
input = str2double(get(hObject,'string'));
handles.PK.sm(handles.i:end)=input;
handles.PK.Filter(handles.i);
handles.PK.CalculateParameters(handles.i);
% handles.PK.Thresh(handles.i);
% handles.PK.Calculate(handles.i);
cla(handles.axes_image)
pos=handles.PK.posper(:,:,handles.i);
M=handles.PK.Mper(:,:,handles.i);
hold(handles.axes_image,'on');
plot(handles.PK.matrix_rough_fluorescences(:,handles.i),'color',handles.col,'parent',handles.axes_image)
plot(handles.PK.matrix_filtered_fluorescences(:,handles.i),'linewidth',2,'color','b','parent',handles.axes_image)
plot(handles.PK.xmc(:,handles.i),handles.PK.mmvg(:,handles.i),'+c','parent',handles.axes_image)
plot(handles.PK.xMc(:,handles.i),handles.PK.M(:,handles.i),'+y','parent',handles.axes_image)
plot(handles.PK.xmr(:,handles.i),handles.PK.mmvd(:,handles.i),'+g','parent',handles.axes_image)
plot(handles.PK.xMr(:,handles.i),handles.PK.M(:,handles.i),'+r','parent',handles.axes_image)
plot(handles.PK.posm(:,handles.i),handles.PK.ms(:,handles.i),'+k','parent',handles.axes_image)
plot(handles.PK.posM(:,handles.i),handles.PK.Ms(:,handles.i),'+k','parent',handles.axes_image)
plot(pos(:),M(:),'+m','parent',handles.axes_image)
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function smooth_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smooth_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function prop_Callback(hObject, eventdata, handles)
% hObject    handle to prop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prop as text
%        str2double(get(hObject,'String')) returns contents of prop as a double
input = str2double(get(hObject,'string'));
handles.PK.prop(handles.i:end)=input;
handles.PK.Filter(handles.i);
handles.PK.CalculateParameters(handles.i);
% handles.PK.Thresh(handles.i);
% handles.PK.Calculate(handles.i);
cla(handles.axes_image)
pos=handles.PK.posper(:,:,handles.i);
M=handles.PK.Mper(:,:,handles.i);
hold(handles.axes_image,'on');
plot(handles.PK.matrix_rough_fluorescences(:,handles.i),'color',handles.col,'parent',handles.axes_image)
plot(handles.PK.matrix_filtered_fluorescences(:,handles.i),'linewidth',2,'color','b','parent',handles.axes_image)
plot(handles.PK.xmc(:,handles.i),handles.PK.mmvg(:,handles.i),'+c','parent',handles.axes_image)
plot(handles.PK.xMc(:,handles.i),handles.PK.M(:,handles.i),'+y','parent',handles.axes_image)
plot(handles.PK.xmr(:,handles.i),handles.PK.mmvd(:,handles.i),'+g','parent',handles.axes_image)
plot(handles.PK.xMr(:,handles.i),handles.PK.M(:,handles.i),'+r','parent',handles.axes_image)
plot(handles.PK.posm(:,handles.i),handles.PK.ms(:,handles.i),'+k','parent',handles.axes_image)
plot(handles.PK.posM(:,handles.i),handles.PK.Ms(:,handles.i),'+k','parent',handles.axes_image)
plot(pos(:),M(:),'+m','parent',handles.axes_image)
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function prop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
