clc
clear
close all;
folder='/data1/thoman/ownCloud/flux_calcique/Comparaison Matlab-Excel P Seyer/';
%% load data
    [file,rough_data_foldername]=uigetfile([folder,'*.*']);
    rough_data_pathname=[rough_data_foldername, file];
 if strcmp(file(end-2:end),'csv')
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
results_foldername=[rough_data_foldername,filesep,'Results_',file(1:end-5), filesep];
if ~exist(results_foldername,'file')
    mkdir(results_foldername);
end

PK=AnalysisPeaks(matrix_rough_data,'param_filter',pol_length,'smoothness',sm,'prop',prop,'type',type,'baselinefit',1,'bool_baselineref',1);


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



