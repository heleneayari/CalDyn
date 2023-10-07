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

% Last Modified by GUIDE v2.5 06-Oct-2023 22:02:51

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

function PeakAnalysis_OpeningFcn(hObject, ~, handles, varargin)
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
handles.bb=0;

set(handles.frame_length,'string',num2str(handles.PK.vector_filtering_frame_length(handles.i)))
set(handles.prop,'string',num2str(handles.PK.prop(handles.i)))
set(handles.type,'value',handles.PK.type)
set(handles.th_smpks,'string',num2str(handles.PK.th_smpks(handles.i)))
set(handles.th_medpks,'string',num2str(handles.PK.th_medpks(handles.i)))
set(handles.fac_multi,'string',num2str(handles.PK.th_multi(handles.i)))
% set(handles.smooth_length,'string',num2str(handles.PK.sm(handles.i)))
if handles.PK.type==2
    
    set(handles.panel_stat,'Visible','off');
    error=1;
%     while error
%         try
            handles.PK.Filter(handles.i);
            handles.PK.CalculateParameters(handles.i);
            error=0;
            set(handles.smooth_length,'string',num2str(handles.PK.sm(handles.i)))
%         catch
%             handles.PK.sm(handles.i)=handles.PK.sm(1)+10;
%             
%         end
%         
%     end
else
    error=1;
    while error
        try
            handles.PK.Filter(handles.i);
            handles.PK.CalculateParameters(handles.i);
            error=0;
            set(handles.smooth_length,'string',num2str(handles.PK.sm(handles.i)))
        catch
            handles.PK.sm(:)=handles.PK.sm(1)+1;
            
        end
    end
    
    error=1;
    while error&&handles.PK.prop(handles.i)>0
        try
            
            handles.PK.CalculateParameters(handles.i);
            handles.PK.prop(handles.i)=handles.PK.prop(handles.i)-0.01;
            set(handles.prop,'string',num2str(handles.PK.prop(handles.i),'%0.3f'))
        catch
            handles.PK.prop(handles.i)=handles.PK.prop(handles.i)+0.01;
            handles.PK.prop(handles.i)
            handles.PK.CalculateParameters(handles.i);
            set(handles.prop,'string',num2str(handles.PK.prop(handles.i),'%0.3f'))
            error=0;
        end
    end
    
    
end
plot_graphs(handles);

guidata(hObject, handles);
uiwait(handles.figure1);

function varargout = PeakAnalysis_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.PK;
if (isfield(handles,'closeFigure') && handles.closeFigure)
    figure1_CloseRequestFcn(hObject, eventdata, handles)
end
function figure1_CloseRequestFcn(hObject, ~, ~)
delete(hObject);

function Discard_Callback(hObject, ~, handles)

guidata(hObject, handles);

function pushbutton_save_Callback(hObject, ~, handles)
hgexport(handles.figure1, [handles.ResDir, 'Peak_Analysis_cell', num2str(handles.i),'.png'], hgexport('factorystyle'), 'Format', 'png');


handles.closeFigure = true;
uiresume(handles.figure1);
guidata(hObject, handles);






function frame_length_Callback(hObject, ~, handles)
% hObject    handle to frame_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_length as text
%        str2double(get(hObject,'String')) returns contents of frame_length as a double
input = str2double(get(hObject,'string'));

old=handles.PK.vector_filtering_frame_length(handles.i);

try
    if input>2 && mod(input,2)==1
        handles.PK.vector_filtering_frame_length(handles.i:end)=input;
    else
        if input<3
            handles.PK.vector_filtering_frame_length(handles.i:end)=3;
            set(handles.frame_length,'string',num2str(3));
        else
            handles.PK.vector_filtering_frame_length(handles.i:end)=input+1;
            set(handles.frame_length,'string',num2str(input+1));
        end
    end
    handles.PK.Filter(handles.i);
    if handles.bb
        handles.PK.remove_base(handles.i);
    end
    handles.PK.CalculateParameters(handles.i);
    
    
    plot_graphs(handles);
    
    
catch
    set(handles.frame_length,'string',num2str(old))
    handles.PK.vector_filtering_frame_length(handles.i:end)=old;
    handles.PK.Filter(handles.i);
end


guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function frame_length_CreateFcn(hObject, ~, ~)
% hObject    handle to frame_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function smooth_length_Callback(hObject, ~, handles)
% hObject    handle to smooth_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of smooth_length as text
%        str2double(get(hObject,'String')) returns contents of smooth_length as a double
input = str2double(get(hObject,'string'));
old=handles.PK.sm(handles.i);
% try
    handles.PK.sm(handles.i)=input;
    handles.PK.Filter(handles.i);
    if handles.bb
        handles.PK.remove_base(handles.i);
    end
    handles.PK.CalculateParameters(handles.i);
    plot_graphs(handles);
    
    
% catch
%     set(handles.smooth_length,'string',num2str(old))
%     handles.PK.sm(handles.i)=old;
%     handles.PK.Filter(handles.i);
%     
%     
% end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function smooth_length_CreateFcn(hObject, ~, ~)
% hObject    handle to smooth_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function prop_Callback(hObject, ~, handles)
% hObject    handle to prop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prop as text
%        str2double(get(hObject,'String')) returns contents of prop as a double
input = str2double(get(hObject,'string'));
old=handles.PK.prop(handles.i);
try
    
    handles.PK.prop(handles.i)=input;
    handles.PK.CalculateParameters(handles.i);
    
    plot_graphs(handles);
    
catch
    handles.PK.prop(handles.i)=old;
    set(handles.prop,'string',old);
    handles.PK.Filter(handles.i);
    
    
end


guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function prop_CreateFcn(hObject, ~, ~)
% hObject    handle to prop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function handles=plot_graphs(handles)
%for plots

pos=handles.PK.posper(:,:,handles.i);
M=handles.PK.Mper(:,:,handles.i);

cla(handles.axes_image)
pos=handles.PK.posper(:,:,handles.i);
M=handles.PK.Mper(:,:,handles.i);
hold(handles.axes_image,'on');

plot(handles.PK.vector_time,handles.PK.matrix_rough_fluorescences(:,handles.i),'color',handles.col,'parent',handles.axes_image)
plot(handles.PK.vector_time,handles.PK.matrix_filtered_fluorescences(:,handles.i),'linewidth',2,'color','b','parent',handles.axes_image)
%plot(handles.PK.vector_time,handles.PK.base(:,handles.i),'linewidth',1,'color','r','parent',handles.axes_image)
plot(handles.PK.xmc(:,handles.i),handles.PK.mmvg(:,handles.i),'+c','linewidth',2,'parent',handles.axes_image)
plot(handles.PK.xMc(:,handles.i),handles.PK.M(:,handles.i),'+y','linewidth',2,'parent',handles.axes_image)
plot(handles.PK.xmr(:,handles.i),handles.PK.mmvd(:,handles.i),'+g','linewidth',2,'parent',handles.axes_image)
plot(handles.PK.xMr(:,handles.i),handles.PK.M(:,handles.i),'+r','linewidth',2,'parent',handles.axes_image)
plot(handles.PK.posm(:,handles.i),handles.PK.ms(:,handles.i),'+k','linewidth',2,'parent',handles.axes_image)
plot(handles.PK.posM(:,handles.i),handles.PK.Ms(:,handles.i),'+k','linewidth',2,'parent',handles.axes_image)
plot(pos(:),M(:),'+m','linewidth',2,'parent',handles.axes_image)
xlabel(handles.axes_image,'Time')

if handles.PK.type==1
    
    handles.PK.statpks(handles.i);
    %     plot(handles.PK.vector_time,(handles.PK.base(:,handles.i)),'color','m','parent',handles.axes_image)
    %     plot(handles.PK.vector_time,(handles.PK.basefit(:,handles.i)),'color','m','parent',handles.axes_image)
    set(handles.Npks_val,'string',num2str(handles.PK.N(handles.i)));
    set(handles.fsmpk_val,'string',num2str(handles.PK.f_smpks(handles.i)));
    set(handles.fmedpk_val,'string',num2str(handles.PK.f_medpks(handles.i)));
    set(handles.fmultipk_val,'string',num2str(handles.PK.f_multipks(handles.i)));

    for ii=1:handles.PK.N(handles.i)

        if handles.PK.ind_medpks(ii,handles.i)
            
            y1=handles.PK.matrix_filtered_fluorescences(round(handles.PK.xmcr(ii,handles.i)):round(handles.PK.posmr(ii,handles.i)),handles.i);            
            x1=handles.PK.vector_time(round(handles.PK.xmcr(ii,handles.i)):round(handles.PK.posmr(ii,handles.i)));
            
            plot(x1,y1,'color',[0.47 0.25 0.80],'linewidth',2,'parent',handles.axes_image)
        end
        
         if handles.PK.ind_smpks(ii,handles.i)
            
            y1=handles.PK.matrix_filtered_fluorescences(round(handles.PK.xmcr(ii,handles.i)):round(handles.PK.posmr(ii,handles.i)),handles.i);            
            x1=handles.PK.vector_time(round(handles.PK.xmcr(ii,handles.i)):round(handles.PK.posmr(ii,handles.i)));
            
            plot(x1,y1,'y','linewidth',2,'parent',handles.axes_image)
         end
        
          if handles.PK.ind_multipks(ii,handles.i)
            y1=handles.PK.matrix_filtered_fluorescences(round(handles.PK.xmcr(ii,handles.i)):round(handles.PK.posmr(ii,handles.i)),handles.i);            
            x1=handles.PK.vector_time(round(handles.PK.xmcr(ii,handles.i)):round(handles.PK.posmr(ii,handles.i)));
            
            plot(x1,y1,'--g','linewidth',2,'parent',handles.axes_image)
         end
        
         
        
    end
    
end

cla(handles.axes2,'reset')

hold(handles.axes2,'on')
yyaxis(handles.axes2,'left')
plot(handles.PK.vector_time,handles.PK.matrix_filtered_fluorescences(:,handles.i),'--','linewidth',1,'color','b','parent',handles.axes2)
plot(handles.PK.vector_time,handles.PK.smooth_signal(:,handles.i),'-','linewidth',1,'color','k','parent',handles.axes2)

yyaxis(handles.axes2,'right')

plot(handles.PK.vector_time,handles.PK.dd2(:,handles.i),'linewidth',2,'color','r','parent',handles.axes2)
plot([0, handles.PK.vector_time(end)],[handles.PK.mder(handles.i), handles.PK.mder(handles.i)],'r','parent',handles.axes2)
plot([0, handles.PK.vector_time(end)],[handles.PK.Mder(handles.i), handles.PK.Mder(handles.i)],'r','parent',handles.axes2)
xlabel('Time')


% --- Executes on selection change in type.
function type_Callback(~, ~, ~)
% hObject    handle to type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from type


% --- Executes during object creation, after setting all properties.
function type_CreateFcn(hObject, ~, ~)
% hObject    handle to type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function th_smpks_Callback(hObject, ~, handles)
% hObject    handle to th_smpks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of th_smpks as text
%        str2double(get(hObject,'String')) returns contents of th_smpks as a double
input = str2double(get(hObject,'string'));
handles.PK.th_smpks(handles.i)=input;
plot_graphs(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function th_smpks_CreateFcn(hObject, ~, ~)
% hObject    handle to th_smpks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function th_medpks_Callback(hObject, ~, handles)
% hObject    handle to th_medpks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of th_medpks as text
%        str2double(get(hObject,'String')) returns contents of th_medpks as a double
input = str2double(get(hObject,'string'));
handles.PK.th_medpks(handles.i)=input;
plot_graphs(handles)
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function th_medpks_CreateFcn(hObject, ~, ~)
% hObject    handle to th_medpks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rem_base_line.
function rem_base_line_Callback(hObject, ~, handles)
% hObject    handle to rem_base_line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.bb=1;
handles.PK.remove_base(handles.i);
handles.PK.CalculateParameters(handles.i);
plot_graphs(handles)
guidata(hObject, handles);




function fac_multi_Callback(hObject, ~, handles)
% hObject    handle to fac_multi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fac_multi as text
%        str2double(get(hObject,'String')) returns contents of fac_multi as a double
input = str2double(get(hObject,'string'));
handles.PK.th_multi(handles.i)=input;
plot_graphs(handles)
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function fac_multi_CreateFcn(hObject, ~, ~)
% hObject    handle to fac_multi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
