clear; close all; clc;
clc
clear
close all;
folder=['Data' filesep 'optical_mapping' filesep];


%% load image
[file,foldername]=uigetfile([folder,'*.*']);
pathname=[foldername, file];


%% We will read the image make an excel file with a few points and then run it
%using the usual code

info=imfinfo(pathname);
T=length(info);

A=imread(pathname,1);
%% we wil make a mask before detecting some random points in the image
pos=strfind(file,'.');
if ~exist([foldername,file(1:pos-1),'_bw.mat'],'file')
    figure;
    title('draw ROI')
    imagesc(A)
    colormap(gray)
    axis equal
    bw=roipoly;
    save([foldername,file(1:pos-1),'_bw.mat'],'bw')
    
else
    load([foldername,file(1:pos-1),'_bw.mat'])
end
%% choose points
clear pts ptst
figure('Name','ptssuperimposed')
hold on
figToolbarFix
imagesc(A)
colormap(gray)
axis equal
axis tight
axis off

[ptst(:,1), ptst(:,2)] = localMaximum_h(double(A), 10, -1, 1000);



ind = sub2ind(size(bw), ptst(:,1), ptst(:,2));
ind2 = bw(ind)>0;
pts(:,:) = ptst(ind2,:);
plot(pts(:,2),pts(:,1),'r+')
set(gca,'Ydir','reverse')
indf = sub2ind(size(bw), pts(:,1), pts(:,2));

%% load the whole tiff
for ii=1:T
    
    imf(:,:,ii)=imread(pathname,ii);
    imff(:,:,ii)=medfilt2(imf(:,:,ii),[3 3]);
end

%%
clear titi

for uu=1:length(pts)
    
    values(:,uu)=[imff(pts(uu,1),pts(uu,2),:)];
    
end
%%
time = (1:T)';


valuesf=cat(2,time,values);
myTable = array2table(valuesf);

%% Add column names
nomsColonnes =cellstr([{'Time'}, strcat('Column_', string(1:length(pts)))]);
myTable.Properties.VariableNames = nomsColonnes;

%% write table to file in excel format
nomFichier =[foldername,file(1:pos-1), '.xlsx'];
writetable(myTable, nomFichier);

disp(['Timeseries saved in file :', nomFichier]);




