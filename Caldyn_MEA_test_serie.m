clc
clear
close all;

% folder='/data1/thoman/Lamia/';
folder='/data1/thoman/Nextcloud_cnrs/Albano/Text Data Experiment 1/';
%folder='D:\Helene\nextcloud_cnrs\Albano\Text Data Experiment 1\'
% folder='C:\Users\HEDY\NextCloud_cnrs\Albano\Text Data Experiment 1\'

%% load data


files=dir([folder,'*textdata.txt']);


%% param 
    sm=500;
    prop=0.8; 
    props=0.1;
    type=2;%2 for MEA;1 for the rest
    param_filter=300;
    col_line=1; %  set to 1, to get each parameter on  a different column in the excels files
    list_param={'N_pks','Amp_mea','Tau_Start_Peak','Tau_Start_End'};



for ff=5:length(files)
    file=files(ff).name;
    rough_data_pathname=[folder, file];

       warning('off','MATLAB:table:ModifiedAndSavedVarnames')

       matrix_rough_data=table2array(readtable(rough_data_pathname));



results_foldername=[folder,filesep,'Results_',file(1:end-4), filesep];
if ~exist(results_foldername,'file')
    mkdir(results_foldername);
end

%%
close all;
clc;

PK=AnalysisPeaks(matrix_rough_data,'prop',prop,'props',props,'type',type,'param_filter',param_filter,'Smoothness',sm,'list_param_name',list_param,'col_line',col_line);





for i=1:PK.number_cells

 PK=PeakAnalysis(PK,i,results_foldername);

end


props=PK.props(i);
%%

    results_pathname=[results_foldername filesep 'Results_test' file(1:end-4),'.xlsx'];
   % save([results_foldername filesep 'Results_' file(1:end-4),'.mat'],'PK','-v7.3');
    PK.Save(results_pathname);
    
end% 
