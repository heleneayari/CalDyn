clc
clear
% close all;
folder='/data1/thoman/ownCloud/Albano/Yvonne/';
folder='/data1/thoman/ownCloud/flux_calcique/Signaux_calciques/'
% folder='/data1/thoman/ownCloud/Albano/wetransfer_excel-files_2024-02-05_1416/Excel files/'
% folder='/data1/thoman/ownCloud/Albano/'
%folder='C:\Users\HEDY\ownCloud\Albano\wetransfer_excel-files_2024-02-05_1416\'
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
        matrix_rough_data=table2array(readtable(rough_data_pathname,'Delimiter',','));
    end
    
else
    
    %      matrix_rough_data=xlsread(rough_data_pathname);
    warning('off','MATLAB:table:ModifiedAndSavedVarnames')
    
    matrix_rough_data=table2array(readtable(rough_data_pathname));
    
    if iscell(matrix_rough_data)
        
        clear matrix_rough_data
        tutu=readtable(rough_data_pathname);
        matrix_rough_data=str2double(strrep(tutu{:,:},',','.'));        
    end
    
end



%% param
pol_length=11;
sm=51;
prop=0.2;
baselinefit=0;%0 for no baseline removal, 1 is for fitting baseline with a constant, 2 with poly1, 3 with poly2, 4 with poly3
type=1;%1 for calcic signals, 2 for electrics
pks_class=0;

th_smpks=0.2;
th_medpks=0.5;
th_multi=0.2;

col_line=1; %  to get each parameter on  a different column in the excels files if set to 1
%%

results_foldername=[rough_data_foldername,filesep,'Results_',file(1:end-5), filesep];
if ~exist(results_foldername,'file')
    mkdir(results_foldername);
end


%%

PK=AnalysisPeaks(matrix_rough_data,'param_filter',pol_length,'smoothness',sm,'prop',prop,'type',type,'pks_class',pks_class,'th_smpks',th_smpks,'th_medpks',th_medpks,'th_multi',th_multi,'baselinefit',baselinefit,'col_line',col_line);
%% test plot
% figure
% plot(PK.vector_time,PK.matrix_rough_fluorescences(:,2))


%% rest



for i=1:PK.number_cells
    
    PK=PeakAnalysis(PK,i,results_foldername);
    close all
    
end

%%



results_pathname=[results_foldername filesep 'Results_' file(1:end-5),'.xlsx'];
save([results_foldername filesep 'Results_' file(1:end-5),'.mat'],'PK');
PK.Save(results_pathname);



