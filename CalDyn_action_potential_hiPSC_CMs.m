clc
clear
close all;
folder=['Data' filesep 'action_potential_hiPSC_CMs' filesep];

%% load data
    [file,foldername]=uigetfile([folder,'*.*']);
    pathname=[foldername, file];

    warning('off','MATLAB:table:ModifiedAndSavedVarnames')
    matrix_rough_data=table2array(readtable(pathname));




%% param 
    pol_length=3;
    sm=50;
    prop=0.05; 


%%
close all
results_foldername=[foldername,filesep,'Results_',file(1:end-5), filesep];
if ~exist(results_foldername,'file')
    mkdir(results_foldername);
end
PK=AnalysisPeaks(matrix_rough_data,'param_filter',pol_length,'smoothness',sm,'prop',prop);

%% rest


% 
for i=1:PK.number_cells

PK=PeakAnalysis(PK,i,results_foldername);
close all

end

%%

    

    results_pathname=[results_foldername filesep 'Results_' file(1:end-5),'.xlsx'];
    save([results_foldername filesep 'Results_' file(1:end-5),'.mat'],'PK');
    PK.Save(results_pathname);

