clc
clear
close all;
folder='/data1/thoman/ownCloud/Lamia/';
folder='/data1/thoman/Lamia/';
folder='/data1/thoman/ownCloud/Albano/Text Data Experiment 1/'

%% load data
tic
    [file,rough_data_foldername]=uigetfile([folder,'*.*']);
    rough_data_pathname=[rough_data_foldername, file];
%  if strcmp(file(end-2:end),'csv')
%      
%     tutu=readtable(rough_data_pathname,'Delimiter',',');
%     for ii=1:size(tutu,2)
%         matrix_rough_data(1:size(tutu,1)-2,ii)=str2double(strrep(tutu{3:end,ii},',','.'));     
%     end
%     
%  else
       warning('off','MATLAB:table:ModifiedAndSavedVarnames')
%        matrix_rough_datat=table2array(readtable(rough_data_pathname));
       matrix_rough_data=table2array(readtable(rough_data_pathname));
%      matrix_rough_data=xlsread(rough_data_pathname);
%  end
 toc
 

%% param 
    sm=1;
    prop=0.5; 
    type=2;%2 for MEA;1 for the rest
    param_filter=300;
    col_line=1; %  set to 1, to get each parameter on  a different column in the excels files
    list_param={'N_pks','Amp_mea','Tau_mea','BI','FPD'};

results_foldername=[rough_data_foldername,filesep,'Results_',file(1:end-4), filesep];
if ~exist(results_foldername,'file')
    mkdir(results_foldername);
end
%% keep only good columns

% % goodcolumn=[1,2:3:size(matrix_rough_datat,2)];
% goodcolumn=[2,3:2:size(matrix_rough_datat,2)];
% matrix_rough_data=matrix_rough_datat(:,goodcolumn);
% matrix_rough_data=matrix_rough_datat(:,1:end-1);
% cut=round(length(matrix_rough_datat(:,1))/100);
% matrix_rough_data=matrix_rough_datat(cut:end,[1,3:10]);
%% test figure to see if data are correct for debugging
% figure
% 
% plot(matrix_rough_data(:,1),matrix_rough_data(:,2))
% pause


%%
close all;
clc;
PK=AnalysisPeaks(matrix_rough_data,'prop',prop,'type',type,'param_filter',param_filter,'Smoothness',sm,'list_param_name',list_param,'col_line',col_line);





for i=1:PK.number_cells

 PK=PeakAnalysis(PK,i,results_foldername);

end

%%

    results_pathname=[resu<lts_foldername filesep 'Results_test' file(1:end-4),'.xlsx'];
   % save([results_foldername filesep 'Results_' file(1:end-4),'.mat'],'PK','-v7.3');
    PK.Save(results_pathname);


% 
