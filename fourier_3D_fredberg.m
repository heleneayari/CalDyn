clear
close all
clc
%% 
dossier='D:\owncloud\Albano\CTL\'
dossier='/data1/thoman/Albano/transfer_6138471_files_2435448d/'
file=dir([dossier,'*.tif']);
i=6;
nom=file(i).name;

savedossier=[dossier,'resultsfourier/'];
 %savedossier='/data1/thoman/Nextcloud/Sarah/resultsfourier';
if~exist(savedossier,'dir')
    mkdir(savedossier)
end
% cd(savedossier)
%nom='Temp_Hcc70 60';

info=imfinfo([dossier,nom]);
T=length(info);
T=100;
DeltaT=5;%min
sc=0.645;%mu/pixels
%T=length(info);
if mod(T,2)==0
omega=[-T/2:1:T/2-1]/T;
else
    omega=[-(T-1)/2:1:(T-1)/2]/T;
end
ideb=100;
for ii=1:T
    %for ii=1:length(info)
    ima=imread([dossier,nom],ii);
%im(:,:,ii)=imcrop(ima,[1 1 754 754]);
im(:,:,ii)=ima;
end
%% Fit des lorentziennes
% Io=3;
% Gao=1;
% Lo=@(Om,Ga,x) Om*Ga^3./((x.^2-Om^2).^2+x.^2*Ga.^2)+Io*(1/2*Gao)./(x.^2+(1/2*Gao)^2);
% Lo=@(Om,Ga,Io,Gao,x) Om*Ga^3./((x.^2-Om^2).^2+x.^2*Ga.^2)+Io*(1/2*Gao)./(x.^2+(1/2*Gao)^2);
 %Lo=@(I,Om,Gamma,Io,c,x) I*Om*(Gamma)^2./((x.^2-Om^2).^2+x.^2*(Gamma).^2)+Io./(x.^2+(1/2*Io)^2)+c;
  Lo=@(I,Om,Gamma,c,Gao,x) I*Om*(Gamma)^2./((x.^2-Om^2).^2+x.^2*(Gamma).^2)+c-I/2*(1/2*Gao)./(x.^2+(1/2*Gao)^2);
%  Lo=@(Om,Io,c,x) -Om./(x.^2+(1/2*Om)^2)+Io./(x.^2+(1/2*Io)^2)+c;
%Lo=@(Om,Ga,Io,x) Om*Ga^3./((x.^2-Om^2).^2+x.^2*Ga.^2)+Io;
Lo2=@(Io,Gao,x) -Io*(1/2*Gao)./(x.^2+(1/2*Gao)^2);
Lob=@(Om,Ga,x)Om*Ga^3./((x.^2-Om^2).^2+x.^2*Ga.^2);
            %FIA = FourierImageAnalysisModelwphase('image',im);
            FIA = FourierImageAnalysisModel('image',im);
            FIA.performFft;
            %FIA=cutCenter(FIA,8);
            FIA.interpolateFft2D;
          % close all 
            I=2;
            I0=1.4
            Ga0=0.3;
            Om=0.18;
            Gamma=0.18;
            ii=55
            c=2.5;
%             figure
%             hold on
%             plot(omega,Lo(I,Om,Gamma,c,omega))
%             plot(omega,Lo2(I0,Ga0,omega))
%             plot(omega,Lo(I,Om,Gamma,c,omega)+Lo2(I0,Ga0,omega))
%                  plot(omega,FIA.Msz(:,ii))
[X,Y]=meshgrid(FIA.qr(1:end/2)/sc,omega(1:end)/DeltaT);
[X2,Y2]=meshgrid(FIA.qr(1:end)/sc,omega(1:end)/DeltaT);
        W=FIA.Msz(51:end,1:end/2);
  %    W=FIA.Msz(11:end,1:end/2);
       Wt=FIA.Msz(1:end,1:end/2);
            figure
            imagesc(log(FIA.Mszdth))
            [~,I]=max(W);
            tutu=I+50;
          %  tutu=I+10;
           
            freq=omega/DeltaT; %en min-1
%             fu=fit(FIA.qr(1:end/2)'/sc,freq(tutu)','poly1');
            x=FIA.qr(1:end/2)'/sc;
%             y1=fu.p2+fu.p1*x;
%             y2=fu.p2-fu.p1*x;
%             f=fit([1:size(W,2)]',tutu','poly1');
            hf1=figure
            hold on
            surf(X2,Y2,log(FIA.Msz),'edgecolor','none')
            %set(gca,'Xlim',[0 0.22])
            %axis off
            xlabel('k ($\mu m^{-1}$)')
            ylabel('f ($min^{-1}$)')
            %saveas(hf1,[dossier,'logbrillouin.png'])
            %plot(1:size(W,2),51+I(:),'+')
%             plot(f)
%             c=colorbar;
%             c.Limits=[1e9, 1e21];
%             
            hf=figure('Position',[-197   121   731   848]);
            hold on
            surf(X,Y,Wt-max(Wt(:)),'edgecolor','none')
%             plot(x,y1,'k')
            %plot(x,y2,'k')
            shading interp
            set(gca,'Xlim',[0 0.1])
            xlabel('k ($\mu m^{-1}$)')
            ylabel('f ($min^{-1}$)')
          %  saveas(hf,[dossier,'brillouin.png'])
            figure
            imagesc(Wt)
            colorbar
            set(gca,'Clim',[0 10])
            
            

            figure
             plot(FIA.Msz(:,2))
figure
hold  on
            for ii=1:1:10
                plot(FIA.Msz(:,ii))
                
            end
            vec=10:2:40;
            map=jet(length(vec));
            figure
hold  on
cp=1;
            for ii=vec
                plot(FIA.Msz(:,ii),'color',map(cp,:))
     cp=cp+1;
            end
            
         %   close all
            clear Io om
           
            figure
            hold  on
            cp=1;
           
            for ii=49:1:80
                plot(omega,FIA.Msz(:,ii))
                ft = fittype( 'Lo1( x, Om,Ga,c )' );
                          %f=fit(omega',FIA.Mszdth(:,ii),ft,'start',[ 0.05 0.01 5 ])
                           
                        %  f=fit(omega',FIA.Mszdth(:,ii),Lo,'start',[ 1 0.05 1  5 0.1])
                         f=fit(omega',FIA.Msz(:,ii),Lo,'start',[2 0.18  0.18 2.5  0.18])
                plot(f)
%                 plot(omega,Lob(f.Om,f.Om*2,omega))
%                 plot(omega,Lo2(f.Io,f.Io,omega))
                om(cp)=f.Om;
                Io(cp)=f.I;
                Gamma(cp)=f.Gamma;
                Gao(cp)=f.Gao;
                cp=cp+1;
           
            end
            figure
            plot(om)
            figure
            plot(Io)
            figure
            plot(Gamma)
             figure
            plot(Gao)
            figure
            plot(om,Io)
%             figure
%             plot(FIA.Msz(201,:))
            
            %% explore FFT
            
%             
%             figure
%             imagesc(FIA.FftEnergy(:,:,78))
            
            view=FourierSpaceTimeViewer(FIA,nom,savedossier);
              view.Msz
  
  view.Phaseth
  view.Fourierspaceonlyav
  FIA.Image=[];
  FIA.FftImage=[];
  FIA.FftEnergy=[];
  FIA.Fsz=[];

  FIA.phase=[];
  FIA.phaserth=[];
 % save([savedossier,'fourier.mat'],'FIA','-v7')
            