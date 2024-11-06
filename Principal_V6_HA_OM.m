
clc
clear
close all;
folder='/data1/thoman/ownCloud/flux_calcique/Signaux_OM/';
folder='/data1/thoman/Albano/Optical Mapping Datasets/GP Unfiltered 3Hz S1 2020-11-24-121415/test/';
% folder='/data1/thoman/Albano/Optical Mapping Datasets/GP Unfiltered Spontaneous 2020-11-24-121305/test/'
% folder='/data1/thoman/Albano/Optical Mapping Datasets/Mouse Unfiltered 5Hz S1 2022-10-17-145641/test/'

%% add basic_functions to path
filePath = matlab.desktop.editor.getActiveFilename;
pos=strfind(filePath,filesep);
pathloc=filePath(1:pos(end-1));
addpath([pathloc,'basic_functions'])
%% load data
    [file,rough_data_foldername]=uigetfile([folder,'*.*']);
    rough_data_pathname=[rough_data_foldername, file];
 if strcmp(file(end-2:end),'csv')
     try
    tutu=readtable(rough_data_pathname,'Delimiter',';');
    for ii=1:size(tutu,2)
        matrix_rough_data(1:size(tutu,1)-2,ii)=str2double(strrep(tutu{3:end,ii},',','.'));     
    end
     catch
         TT=readtable(rough_data_pathname,'Delimiter',',');
         tab=[1,2:2:size(TT,2)];
         matrix_rough_data=TT{:,tab};
     end
 else
     warning('off','MATLAB:table:ModifiedAndSavedVarnames')
     matrix_rough_data=table2array(readtable(rough_data_pathname));
 end


%% for inverting the matrix 
M=max(matrix_rough_data);
matrix_rough_data(:,2:end)=M(2:end)-matrix_rough_data(:,2:end);
%% set timescale framerate 500Hz;
 matrix_rough_data(:,1)= matrix_rough_data(:,1)/500;
 
 
 %% param 
    pol_length=3;
    sm=30;
    prop=0.3; 
    type=1;%1 for calcic signals, 2 for electrics
    col_line=1; %  set to 1, to get each parameter on  a different column in the excels files
    baselinefit=0; % 0 no fit, 1 fit by a constant, 2 fit by poly1, 3 fit by poly2
    baselineref=0;%1 take fit for value of the baseline 
    auto=0; %1 for treating automatically without any manual intervention all the columns in the file
%%
close all
results_foldername=[rough_data_foldername,filesep,'Results_',file(1:end-5), filesep];
if ~exist(results_foldername,'file')
    mkdir(results_foldername);
end
% PK=AnalysisPeaks(matrix_rough_data,'param_filter',pol_length,'smoothness',sm,'prop',prop,'type',...
%                         type,'col_line',col_line,'baselinefit',baselinefit,'auto',auto);

PK=AnalysisPeaks(matrix_rough_data,'param_filter',pol_length,'smoothness',sm,'prop',prop,...
    'type',type,'col_line',col_line,'baselinefit',baselinefit,'bool_baselineref',baselineref,'auto',auto);

%%
for i=1:PK.number_cells


PK=PeakAnalysis(PK,i,results_foldername);
close all

end

%%

    
% 
    results_pathname=[results_foldername filesep 'Results_' file(1:end-5),'.xlsx'];
    save([results_foldername filesep 'Results_' file(1:end-5),'.mat'],'PK');
    PK.Save(results_pathname);

% %     results_pathname=[results_foldername filesep 'Results_fit_' file(1:end-5),'.xlsx'];
% %     save([results_foldername filesep 'Results_' file(1:end-5),'.mat'],'PK');
% %     PK.Save(results_pathname);

