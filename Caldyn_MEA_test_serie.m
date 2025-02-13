clc
clear
close all;
folder='/data1/thoman/ownCloud/Lamia/';
folder='/data1/thoman/Lamia/';
folder='/data1/thoman/ownCloud/Albano/Text Data Experiment 1/';

% folder='C:\Users\HEDY\ownCloud\Albano\Text Data Experiment 1\';

%% load data
tic

files=dir([folder,'*textdata.txt']);


%% param 
    sm=10;
    prop=0.5; 
    props=0.1;
    type=2;%2 for MEA;1 for the rest
    param_filter=300;
    col_line=1; %  set to 1, to get each parameter on  a different column in the excels files
    list_param={'N_pks','Amp_mea','Tau_mea','BI','FPD'};



for ff=3:length(files)
    file=files(ff).name;
    rough_data_pathname=[folder, file];

       warning('off','MATLAB:table:ModifiedAndSavedVarnames')
%        matrix_rough_datat=table2array(readtable(rough_data_pathname));
       matrix_rough_data=table2array(readtable(rough_data_pathname));

 toc
 



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
