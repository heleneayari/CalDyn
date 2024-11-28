clc
clear
close all;

folder=['Data' filesep 'paced_primary_human_skeletal_myotubes' filesep];


%% load data
    [file,foldername]=uigetfile([folder,'*.*']);
    pathname=[foldername, file];
    
    warning('off','MATLAB:table:ModifiedAndSavedVarnames')
    matrix_rough_data=table2array(readtable(pathname));

    

%% param 
    pol_length=3;
    sm=1;
    prop=0.3; 
    baselinefit=2;% 1 is for fitting baseline with a constant, 2 with poly1, 3 with poly2, 4 with poly3
%%

results_foldername=[foldername,filesep,'Results_',file(1:end-5), filesep];
if ~exist(results_foldername,'file')
    mkdir(results_foldername);
end

PK=AnalysisPeaks(matrix_rough_data,'param_filter',pol_length,'smoothness',sm,'prop',prop,'baselinefit',baselinefit,'bool_baselineref',1);


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



