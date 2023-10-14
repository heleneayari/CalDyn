clc
clear
close all;
folder='/data1/thoman/ownCloud/flux_calcique/';
%% load data
    [file,rough_data_foldername]=uigetfile([folder,'*.*']);
    rough_data_pathname=[rough_data_foldername, file];
 if strcmp(file(end-2:end),'csv')
      'coucou'
    tutu=readtable(rough_data_pathname,'Delimiter',',');
  
    for ii=1:size(tutu,2)
        matrix_rough_data(1:size(tutu,1)-2,ii)=str2double(strrep(tutu{3:end,ii},',','.'));     
    end
 else
%      matrix_rough_data=xlsread(rough_data_pathname);
     matrix_rough_data=table2array(readtable(rough_data_pathname));
 end
% pas=100
% tab=1:pas:length(matrix_rough_data);
% figure
% plot(matrix_rough_data(tab,1),matrix_rough_data(tab,2))
    
rough_data_pathname
%% param 
    pol_length=2001;
    sm=1;
    prop=0.05; 
    type=2;%1 for calcic signals, 2 for electrics ; 3 MEA
    cut_freq=30;

results_foldername=[rough_data_foldername,filesep,'Results_',file(1:end-5), filesep];
if ~exist(results_foldername,'file')
    mkdir(results_foldername);
end

matrix_rough_data(:,2:end)=-matrix_rough_data(:,2:end);
PK=AnalysisPeaks(matrix_rough_data,'Pol_length',pol_length,'smoothness',sm,'prop',prop,'type',type,'cut_freq',cut_freq);



% figure
% plot(PK.vector_time(:,1),PK.matrix_rough_fluorescences(:,1))




for i=1:PK.number_cells
%     for i=3


 PK=PeakAnalysis(PK,i,results_foldername);
% close all

end

%%

    

    results_pathname=[results_foldername filesep 'Results_' file(1:end-5),'.xlsx'];
    save([results_foldername filesep 'Results_' file(1:end-5),'.mat'],'PK');
    PK.Save(results_pathname);


% 
