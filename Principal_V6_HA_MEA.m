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
     matrix_rough_data=table2array(readtable(rough_data_pathname));
 end

%% param 
    sm=1;
    prop=0.4; 
    type=3;%1 for calcic signals, 2 for electrics ; 3 MEA
    param_filter=30;

results_foldername=[rough_data_foldername,filesep,'Results_',file(1:end-5), filesep];
if ~exist(results_foldername,'file')
    mkdir(results_foldername);
end

matrix_rough_data(:,2:end)=-matrix_rough_data(:,2:end);
PK=AnalysisPeaks(matrix_rough_data,'prop',prop,'type',type,'param_filter',param_filter,'Smoothness',sm);





for i=1:PK.number_cells
    
 PK=PeakAnalysis(PK,i,results_foldername);

end

%%

    

    results_pathname=[results_foldername filesep 'Results_' file(1:end-5),'.xlsx'];
    save([results_foldername filesep 'Results_' file(1:end-5),'.mat'],'PK');
    PK.Save(results_pathname);


% 
