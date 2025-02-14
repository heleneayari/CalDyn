clc
clear
close all;
folder=['Data' filesep 'optical_mapping' filesep];


%% load data
    [file,foldername]=uigetfile([folder,'*.*']);
    pathname=[foldername, file];

     warning('off','MATLAB:table:ModifiedAndSavedVarnames')
     matrix_rough_data=table2array(readtable(pathname));




%% set timescale framerate 500Hz;
 matrix_rough_data(:,1)= matrix_rough_data(:,1)/500;
 
 
 %% param 
    pol_length=3;
    sm=30;
    prop=0.3; 

    col_line=1; %  set to 1, to get each parameter on  a different column in the excels files
    baselinefit=0; % 0 no fit, 1 fit by a constant, 2 fit by poly1, 3 fit by poly2
    baselineref=0;%1 take fit for value of the baseline 
    auto=0; %1 for treating automatically without any manual intervention all the columns in the file
%%
close all
results_foldername=[foldername,filesep,'Results_',file(1:end-5), filesep];
if ~exist(results_foldername,'file')
    mkdir(results_foldername);
end
% PK=AnalysisPeaks(matrix_rough_data,'param_filter',pol_length,'smoothness',sm,'prop',prop,'type',...
%                         type,'col_line',col_line,'baselinefit',baselinefit,'auto',auto);

PK=AnalysisPeaks(matrix_rough_data,'param_filter',pol_length,'smoothness',sm,'prop',prop,...
    'col_line',col_line,'baselinefit',baselinefit,'bool_baselineref',baselineref,'auto',auto);

%%
for i=1:PK.number_cells


PK=PeakAnalysis(PK,i,results_foldername);
close all

end

%%

    
% 
    results_pathname=[results_foldername filesep 'Results_' file(1:end-5),'.xlsx'];
    save([results_foldername filesep 'Results_' file(1:end-5),'.mat'],'PK');
    PK.Save(results_pathname);

% %     results_pathname=[results_foldername filesep 'Results_fit_' file(1:end-5),'.xlsx'];
% %     save([results_foldername filesep 'Results_' file(1:end-5),'.mat'],'PK');
% %     PK.Save(results_pathname);
