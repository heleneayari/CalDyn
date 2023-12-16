clc
clear
close all;
folder='/data1/thoman/ownCloud/flux_calcique/Signaux_MED/';

%% load data
    [file,rough_data_foldername]=uigetfile([folder,'*.*']);
    rough_data_pathname=[rough_data_foldername, file];
 if strcmp(file(end-2:end),'csv')

    tutu=readtable(rough_data_pathname,'Delimiter',',');
  
    for ii=1:size(tutu,2)
        matrix_rough_data(1:size(tutu,1)-2,ii)=str2double(strrep(tutu{3:end,ii},',','.'));     
    end
 else
       warning('off','MATLAB:table:ModifiedAndSavedVarnames')
     matrix_rough_datat=table2array(readtable(rough_data_pathname));
%      matrix_rough_data=xlsread(rough_data_pathname);
 end

%% param 
    sm=1;
    prop=0.4; 
    type=2;%2 for MEA;1 for the rest
    param_filter=100;
    list_param={'N_pks','FP_duration','FP_Amp'};

results_foldername=[rough_data_foldername,filesep,'Results_',file(1:end-5), filesep];
if ~exist(results_foldername,'file')
    mkdir(results_foldername);
end
%% keep only good columns

goodcolumn=[1,2:3:size(matrix_rough_datat,2)];
matrix_rough_data=-matrix_rough_datat(:,goodcolumn);

%% test figure to see if data are correct for debugging
figure
plot(matrix_rough_data(:,1),matrix_rough_data(:,2))



%%

PK=AnalysisPeaks(matrix_rough_data,'prop',prop,'type',type,'param_filter',param_filter,'Smoothness',sm,'list_param_name',list_param);





for i=1:PK.number_cells
    
 PK=PeakAnalysis(PK,i,results_foldername);

end

%%

    

    results_pathname=[results_foldername filesep 'Results_' file(1:end-5),'.xlsx'];
    save([results_foldername filesep 'Results_' file(1:end-5),'.mat'],'PK');
    PK.Save(results_pathname);


% 
