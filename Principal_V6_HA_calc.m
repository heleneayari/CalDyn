
clc
clear
% close all;
folder='/data1/thoman/ownCloud/flux_calcique/';
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
     matrix_rough_data=xlsread(rough_data_pathname);
     
 end

    

%% param 
    pol_length=3;
    sm=1;
    prop=0.05; 
    type=1;%1 for calcic signals, 2 for electrics 
    th_smpks=0.2;
    th_medpks=0.5;
    th_multi=0.2;
%%

results_foldername=[rough_data_foldername,filesep,'Results_',file(1:end-5), filesep];
if ~exist(results_foldername,'file')
    mkdir(results_foldername);
end

PK=AnalysisPeaks(matrix_rough_data,'Pol_length',pol_length,'smoothness',sm,'prop',prop,'type',type,'th_smpks',th_smpks,'th_medpks',th_medpks,'th_multi',th_multi);

%% rest



for i=1:PK.number_cells

PK=PeakAnalysis(PK,i,results_foldername);
close all

end

%%

    

    results_pathname=[results_foldername filesep 'Results_' file(1:end-5),'.xlsx'];
    save([results_foldername filesep 'Results_' file(1:end-5),'.mat'],'PK');
    PK.Save(results_pathname);



