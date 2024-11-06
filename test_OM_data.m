clear; close all; clc;
folder='/data1/thoman/Albano/Optical Mapping Datasets/GP Unfiltered 3Hz S1 2020-11-24-121415/';
folder='/data1/thoman/Albano/Optical Mapping Datasets/GP Unfiltered Spontaneous 2020-11-24-121305/'
folder='/data1/thoman/Albano/Optical Mapping Datasets/Mouse Unfiltered 5Hz S1 2022-10-17-145641/'
folders=[folder,'test/'];
if ~exist(folders,'dir')
    mkdir(folders)
end
nom='N256 (IF1-CAM2).tif';

%% We will read the image make an excel file with a few points and then run it
%using the usual code

info=imfinfo([folder,nom]);
T=length(info);

A=imread([folder,nom],1);
% figure;
% imagesc(A)
% colormap(gray)
% axis equal
% bw=roipoly;
% 
% save([folder,'bw.mat'],'bw')
% pause

load([folder,'bw.mat'])

%% choose points
clear pts ptst
figure
hold on
imagesc(bw)
[ptst(:,1), ptst(:,2)] = localMaximum_h(double(A), 10, -1, 1000);



    ind = sub2ind(size(bw), ptst(:,1), ptst(:,2));
    ind2 = bw(ind)>0;
    pts(:,:) = ptst(ind2,:);
plot(pts(:,2),pts(:,1),'r+')    

   indf = sub2ind(size(bw), pts(:,1), pts(:,2));


for ii=1:T
    imf(:,:,ii)=imread([folder,nom],ii);
    imff(:,:,ii)=medfilt2(imf(:,:,ii),[3 3]);
end

%%
clear titi

for uu=1:length(pts)
    
    titi(:,uu)=[imff(pts(uu,1),pts(uu,2),:)];

end
%%
temps = (1:T)';
%%Créer la table
% Combiner les données de temps et les données aléatoires dans une table

tutu=cat(2,temps,titi);
maTable = array2table(tutu);

% Ajouter des noms de colonnes pour chaque colonne de données
nomsColonnes =cellstr([{'Temps'}, strcat('Colonne_', string(1:length(pts)))]);
maTable.Properties.VariableNames = nomsColonnes;

% Étape 3 : Enregistrer la table au format Excel
nomFichier = 'mes_donnees_cam2.xlsx';
writetable(maTable, [folders,nomFichier]);

disp(['La table a été enregistrée dans le fichier ', nomFichier]);




