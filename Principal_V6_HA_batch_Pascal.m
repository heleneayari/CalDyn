clc
clear
close all;
folder='/data1/thoman/ownCloud/flux_calcique/Comparaison Matlab-Excel P Seyer/';
folder='/data1/thoman/ownCloud/pascal/'

%% load  all excell files data

files=[dir([folder,'*.xlsx']);dir([folder,'*.csv'])];

%%
for ff=1:length(files)
    nom=files(ff).name;
    rough_data_pathname=[folder, nom];
 if strcmp(nom(end-2:end),'csv')
    tutu=readtable(rough_data_pathname,'Delimiter',';');
    for ii=1:size(tutu,2)
        matrix_rough_data(1:size(tutu,1)-2,ii)=str2double(strrep(tutu{3:end,ii},',','.'));     
    end
 else
     warning('off','MATLAB:table:ModifiedAndSavedVarnames')
     matrix_rough_data=table2array(readtable(rough_data_pathname));
 end

    

%% param 
    pol_length=3;
    sm=1;
    prop=0.3; 
    type=1;% 2 for MEA, 1 for the rest 
%%
close all
results_foldername=[folder,filesep,'Results_',nom(1:end-5), filesep];
if ~exist(results_foldername,'file')
    mkdir(results_foldername);
end

if exist('param','var')
 PK=AnalysisPeaks(matrix_rough_data,'param_filter',param.param_filter,'smoothness',param.smoothness,'prop',param.prop,'type',param.type,'baselinefit',param.baselinefit,'bool_baselineref',param.bool_baselineref,...
     'list_calc',param.list_calc,'list_param_name',param.list_param_name,'pks_class',param.pks_class);

else
PK=AnalysisPeaks(matrix_rough_data,'param_filter',pol_length,'smoothness',sm,'prop',prop,'type',type,'baselinefit',1,'bool_baselineref',1);
end

%% rest


% 
for i=1:PK.number_cells


PK=PeakAnalysis(PK,i,results_foldername);
close all

end

%%

    

    results_pathname=[results_foldername filesep 'Results_' nom(1:end-5),'.xlsx'];
    save([results_foldername filesep 'Results_' nom(1:end-5),'.mat'],'PK');
    PK.Save(results_pathname);
    param=PK.getallinputparams;

    clear PK

end
