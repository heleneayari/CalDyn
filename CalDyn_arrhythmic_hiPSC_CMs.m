clc
clear
% close all;
% path to example, please replace with your own path to data
folder=['Data' filesep 'arrhythmic_hiPSC_CMs' filesep];

%% load data
[file,foldername]=uigetfile([folder,'*.*']);
pathname=[foldername, file];

warning('off','MATLAB:table:ModifiedAndSavedVarnames')
matrix_rough_data=table2array(readtable(pathname));

% if data are saved with a comma for decimal instead of a point as used in
% matlab
    if iscell(matrix_rough_data)
        
        clear matrix_rough_data
        data_temp=readtable(pathname);
        matrix_rough_data=str2double(strrep(data_temp{:,:},',','.'));        
    end


%% param
pol_length=11;
sm=51;
prop=0.1;
baselinefit=0;%0 for no baseline removal, 1 is for fitting baseline with a constant, 2 with poly1, 3 with poly2, 4 with poly3

pks_class=1;

th_smpks=0.25;
th_medpks=0.5;
th_multi=0.15;

col_line=1; %  to get each parameter on  a different column in the excels files if set to 1
%%

results_foldername=[foldername,filesep,'Results_',file(1:end-5), filesep];
if ~exist(results_foldername,'file')
    mkdir(results_foldername);
end


%%

PK=AnalysisPeaks(matrix_rough_data,'param_filter',pol_length,'smoothness',sm,'prop',prop,'pks_class',pks_class,'th_smpks',th_smpks,'th_medpks',th_medpks,'th_multi',th_multi,'baselinefit',baselinefit,'col_line',col_line);




for i=1:PK.number_cells
    
    PK=PeakAnalysis(PK,i,results_foldername);
    close all
    
end

%%



results_pathname=[results_foldername filesep 'Results_' file(1:end-5),'.xlsx'];
save([results_foldername filesep 'Results_' file(1:end-5),'.mat'],'PK');
PK.Save(results_pathname);



