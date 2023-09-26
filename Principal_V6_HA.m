
clc
clear
folder='/data1/thoman/ownCloud/flux_calcique/';
    %% load data
    [file,rough_data_foldername]=uigetfile([folder,'*.xls*']);
    rough_data_pathname=[rough_data_foldername, file];
    matrix_rough_data=xlsread(rough_data_pathname);
    
slow=0;
if slow
    pol_length=51;
    sm=200;
    prop=0.2;
else
    pol_length=11;
    sm=20;
    prop=0.1;
end
    

%%
close all
results_foldername=[rough_data_foldername,filesep,'Results',filesep];
if ~exist(results_foldername,'file')
    mkdir(results_foldername);
end

PK=AnalysisPeaks(matrix_rough_data,'Pol_length',pol_length,'smoothness',sm,'prop',prop);

%% rest


for i=1:PK.number_cells
    
PK=PeakAnalysis(PK,i,results_foldername);
close all

end

%%

    

    results_pathname=[results_foldername filesep 'Results.xlsx'];
    PK.Save(results_pathname);


close all;
