classdef AnalysisPeaks < handle
    properties
        
        number_acquisitions
        number_cells
        vector_time
        matrix_rough_fluorescences
        vector_filtering_polynomial_order
        vector_filtering_frame_length
        matrix_filtered_fluorescences
        sm
        prop
        PixelSize
        framerate
        SamplingFrequency
        Tauc
        Taur
        Taud
        Aire
        Taus
        pc
        pr
        hr
        hc
        M
        mmvg
        mmvd
        xmc
        xMc
        xmr
        xMr
        locc
        locr
        posm
        posM
        ms
        Ms
        posper
        Mper
        ltab
        dd2
        smooth_signal
 
        
    end
    methods
        function PK = AnalysisPeaks(varargin)
            p = inputParser;
            addRequired(p, 'Signal');
            addParameter(p, 'pkhgt',50)
            addOptional(p, 'PixelSize', 1);
            addOptional(p, 'SamplingFrequency', 1);
            addOptional(p, 'Pol_length', 111);
            addOptional(p, 'Smoothness', 200);
            addOptional(p, 'proportion', 0.1);
            parse(p, varargin{:});

            pol_length=p.Results.Pol_length;
            sm=p.Results.Smoothness;
            PK.ltab=10000;
            prop=p.Results.proportion;
            PK.number_acquisitions=length(p.Results.Signal(~isnan(p.Results.Signal(:,1)),1));
            PK.number_cells=size(p.Results.Signal,2)-1;
            PK.vector_time=p.Results.Signal(1:PK.number_acquisitions,1);
            PK.PixelSize = p.Results.PixelSize;
            PK.framerate=1/median(diff(PK.vector_time));
%             PK.SamplingFrequency=p.Results.SamplingFrequency;
%             PK.framerate = (PK.SamplingFrequency)/1000; %fps in ms


            PK.matrix_rough_fluorescences=p.Results.Signal(1:PK.number_acquisitions,2:(PK.number_cells+1));
            PK.vector_filtering_polynomial_order=2*ones(1,PK.number_cells);
            PK.vector_filtering_frame_length=pol_length*ones(1,PK.number_cells);
            PK.matrix_filtered_fluorescences=zeros(PK.number_acquisitions,PK.number_cells);
            PK.smooth_signal=zeros(PK.number_acquisitions,PK.number_cells);
            PK.dd2=zeros(PK.number_acquisitions,PK.number_cells);

            
            PK.sm=sm*ones(1,PK.number_cells);
            PK.prop=prop*ones(1,PK.number_cells);
            PK.xmr=nan*ones(PK.ltab,PK.number_cells);
            PK.xMr=nan*ones(PK.ltab,PK.number_cells);
            PK.xmc=nan*ones(PK.ltab,PK.number_cells);
            PK.xMc=nan*ones(PK.ltab,PK.number_cells);
            PK.mmvg=nan*ones(PK.ltab,PK.number_cells);
            PK.mmvd=nan*ones(PK.ltab,PK.number_cells);
            PK.M=nan*ones(PK.ltab,PK.number_cells);
            PK.Tauc=nan*ones(PK.ltab,PK.number_cells);
            PK.Taur=nan*ones(PK.ltab,PK.number_cells);
            PK.Taud=nan*ones(PK.ltab,PK.number_cells);
            PK.Taus=nan*ones(PK.ltab,PK.number_cells);
            PK.Tauc=nan*ones(PK.ltab,PK.number_cells);
            PK.Tauc=nan*ones(PK.ltab,PK.number_cells);
            PK.pc=nan*ones(PK.ltab,PK.number_cells);
            PK.pr=nan*ones(PK.ltab,PK.number_cells);
            PK.hc=nan*ones(PK.ltab,PK.number_cells);
            PK.hr=nan*ones(PK.ltab,PK.number_cells);
            PK.Aire=nan*ones(PK.ltab,PK.number_cells);
%             PK.locc=nan*ones(PK.ltab,PK.number_cells);
%             PK.locr=nan*ones(PK.ltab,PK.number_cells);
            PK.posm=nan*ones(PK.ltab,PK.number_cells);
            PK.posM=nan*ones(PK.ltab,PK.number_cells);
            PK.ms=nan*ones(PK.ltab,PK.number_cells);
            PK.Ms=nan*ones(PK.ltab,PK.number_cells);
            PK.posper=nan*ones(PK.ltab,5,PK.number_cells);
            PK.Mper=nan*ones(PK.ltab,5,PK.number_cells);     
            
            
        end
        
        
        function PK=Filter(PK,varargin)
            i=varargin{1};
            PK.matrix_filtered_fluorescences(:,i)=sgolayfilt(PK.matrix_rough_fluorescences(:,i),PK.vector_filtering_polynomial_order(i),PK.vector_filtering_frame_length(i)); % in AU
          %  PK.matrix_filtered_fluorescences(:,i)=PK.matrix_rough_fluorescences(:,i);
        end
        

  
        
        function PK = CalculateParameters(PK,varargin)
            p=inputParser;
            addRequired(p,'i')
            parse(p, varargin{1:end});
            i=p.Results.i;
            dper=[0.8 0.7 0.5 0.3 0.1];
            PK.posm(:,i)=nan*ones(PK.ltab,1);
            PK.posM(:,i)=nan*ones(PK.ltab,1);
            PK.ms(:,i)=nan*ones(PK.ltab,1);
            PK.Ms(:,i)=nan*ones(PK.ltab,1);
            PK.xMr(:,i)=nan*ones(PK.ltab,1);
            PK.xmc(:,i)=nan*ones(PK.ltab,1);
            PK.xMc(:,i)=nan*ones(PK.ltab,1);
            PK.mmvg(:,i)=nan*ones(PK.ltab,1);
            PK.mmvd(:,i)=nan*ones(PK.ltab,1);
            PK.M(:,i)=nan*ones(PK.ltab,1);
            PK.Tauc(:,i)=nan*ones(PK.ltab,1);
            PK.Taur(:,i)=nan*ones(PK.ltab,1);
            PK.Taud(:,i)=nan*ones(PK.ltab,1);
            PK.Taus(:,i)=nan*ones(PK.ltab,1);
            PK.Tauc(:,i)=nan*ones(PK.ltab,1);
            PK.pc(:,i)=nan*ones(PK.ltab,1);
            PK.pr(:,i)=nan*ones(PK.ltab,1);
            PK.hc(:,i)=nan*ones(PK.ltab,1);
            PK.hr(:,i)=nan*ones(PK.ltab,1);
            PK.Aire(:,i)=nan*ones(PK.ltab,1);
%             PK.locc(:,i)=nan*ones(PK.ltab,1);
%             PK.locr(:,i)=nan*ones(PK.ltab,1);

            PK.Mper(:,:,i)=nan*ones(PK.ltab,5,1);
            PK.posper(:,:,i)=nan*ones(PK.ltab,5,1);
          
            
            
           
            
            
            win=3;
            sl=PK.sm(i);
            Signal=PK.matrix_filtered_fluorescences(~isnan(PK.matrix_filtered_fluorescences(:,i)),i);
            %             Signal=    PK.matrix_rough_fluorescences(:,i);
            
            
            ss=smooth(Signal,sl,'loess');
            
            
            
            dd2=gradient(ss);
            PK.dd2(1:length(Signal),i)=dd2;
            PK.smooth_signal(1:length(Signal),i)=ss;
            dd=gradient(Signal);
            
            %                 if(sl~=80||prop~=0.05)
            %                                 figure
            %                                 hold on
            %                            %    plot(dd)
            %                                 plot(dd2)
            %
            %                 end
            mm=max(dd2);
            nn=min(dd2);
            
            [~, locct] = findpeaks(dd2,'MinPeakHeight' ,PK.prop(i)*mm);
            [~, locrt] = findpeaks(-dd2,'MinPeakHeight' ,-PK.prop(i)*nn);
            
            clear maxc locc locr minr
            for ii=1:length(locct)
                [maxc(ii),locc(ii)]=max(dd(max(1,locct(ii)-win):min(locct(ii)+win,length(dd))));
                locc(ii)=max(1,locct(ii)-win)+locc(ii)-1;
            end
            for ii=1:length(locrt)
                [minr(ii),locr(ii)]=min(dd(max(1,locrt(ii)-win):min(locrt(ii)+win,length(dd))));
                locr(ii)=locr(ii)+max(1,locrt(ii)-win)-1;
            end
            
            %                                                                  figure
            %                                                 hold on
            %                                                 plot(dd2)
            %                                                 plot(locc,maxc,'+r')
            %                                                 plot(locr,minr,'+g')
            %                                                 plot(sl,dd2(sl),'+y')
            %                                                 pause
            %
            
            [locr,ia,~]=unique(locr);
            minr=minr(ia);
            [locc,ia,~]=unique(locc);
            maxc=maxc(ia);
            
            
            %% for keeping the first peak in good cases
            
            
            [counts,ed]=histcounts(abs(dd2));
            T=otsuthresh(counts);
            valt=T*max(ed);
            mini=std(abs(dd2(abs(dd2)<valt)));
            
            sl2=1;
            if abs(dd2(sl))<2*mini
                locr=[sl2 locr];
                minr=[dd2(sl2) minr];
                
                
            end
            
            
            
            
            
            %% for keeping last peak in good cases
            if abs(dd2(end-sl))<2*mini
                locc=[locc length(dd2)-sl2];
                maxc=[maxc dd2(end-sl2)];
            end
            
            
            
            
            %
            
            %% calculate parameters for complete peaks
%             figure;
%             hold on
%             plot(Signal,'linewidth',2)
%             plot(PK.matrix_rough_fluorescences(:,i),'color',[0.9 0.9 0.9])
            count=1;
            
            for ii=1:length(locc)-1
                
                pl= find(diff(locr>locc(ii))==1)+1;
                %                     if isempty(pl)
                %                         pl=1;
                %                     end
                if locr(pl)<locc(ii+1)
                    
                    if (pl-1)>0
                        
                        if locc(ii)>locr(pl-1)
               
                            %                                 
                            [M(count),ll(ii)]=max(Signal(locc(ii):locr(pl)));
                            ll(ii)=ll(ii)+locc(ii)-1;
                            
                            [mmvg(count),vvg(ii)]=min(ss(locr(pl-1):locc(ii)));
%                              [mmvg(count),vvg(ii)]=min(Signal(locr(pl-1):locc(ii)));
                            vvg(ii)=vvg(ii)+locr(pl-1)-1;
                            [mmvd(count),vvd(ii)]=min(ss(locr(pl):locc(ii+1)));
%                             [mmvd(count),vvd(ii)]=min(Signal(locr(pl):locc(ii+1)));
                            vvd(ii)=vvd(ii)+locr(pl)-1;
                                                           
                            
                            
                            b=Signal(locc(ii))-maxc(ii)*locc(ii);
                            xmc(count)=(mmvg(count)-b)/maxc(ii);
                            xMc(count)=(M(count)-b)/maxc(ii);
                            b=Signal(locr(pl))-minr(pl)*locr(pl);
                            xMr(count)=(M(count)-b)/minr(pl);
                            xmr(count)=(mmvd(count)-b)/minr(pl);
                            
                            %                                 %
%                             plot(xmc(count),mmvg(count),'+b')
%                             plot(xMc(count),M(count),'+y')
%                             plot(xmr(count),mmvd(count),'+g')
%                             plot(xMr(count),M(count),'+r')
                            
                            PK.Tauc(count,i)=(xMc(count)-xmc(count))/PK.framerate;
                            PK.Taur(count,i)=(xmr(count)-xMr(count))/PK.framerate;
                            PK.Taud(count,i)=(xmr(count)-xMc(count))/PK.framerate;
                            PK.Aire(count,i)=sum(Signal(max(1,round(xmc(count))):min(round(xmr(count)),length(Signal))));
                            PK.Taus(count,i)=(xmr(count)-xmc(count))/PK.framerate;

                 
                            
                            PK.pc(count,i)=maxc(ii)*PK.PixelSize*PK.framerate;
                            PK.pr(count,i)=minr(pl)*PK.PixelSize*PK.framerate;
                            PK.hr(count,i)=-minr(pl)*(xmr(count)-xMr(count))*PK.PixelSize;
                            PK.hc(count,i)=maxc(ii)*(xMc(count)-xmc(count))*PK.PixelSize;
                            PK.M(count,i)=M(count)*PK.PixelSize;
                            PK.mmvg(count,i)=mmvg(count)*PK.PixelSize;
                            PK.mmvd(count,i)=mmvd(count)*PK.PixelSize;

                            [Ms(count),posM(count)]=max(Signal(round(max(1,xmc(count))):round(min(length(Signal),xmr(count)))));
                            posM(count)=posM(count)+xmc(count)-1;

                            
                            count=count+1;
                        end
                    end
                end
            end
            
            for kk=1:length(xMr)-1
                [ms(kk),posm(kk)]=min(Signal(round(xMr(kk)):round(xMr(kk+1))));
                posm(kk)=posm(kk)+xMr(kk)-1;
               
            end
            [ms(length(xMr)),posm(length(xMr))]=min(Signal(round(xMr(length(xMr))):end));
            posm(length(xMr))=posm(length(xMr))+xMr(length(xMr))-1;
         
            for uu=1:5
            for kk=1:length(xMr(:))
                dd=dper(uu)*(Ms(kk)-ms(kk));
                tutu=find(diff(Signal(round(posM(kk)):round(posm(kk)))-ms(kk)>dd)==-1);
                if ~isempty(tutu)
                posper(kk,uu)=tutu(1);
                posper(kk,uu)=posper(kk,uu)+posM(kk)-1;
                Mper(kk,uu)=ms(kk)+dd;
               
                end
            end
            end
            
            PK.posper(1:size(posper,1),:,i)=posper/PK.framerate;
            PK.Mper(1:size(posper,1),:,i)=Mper*PK.PixelSize;
            PK.posm(1:length(xmc),i)=posm/PK.framerate;
            PK.posM(1:length(xmc),i)=posM/PK.framerate;
            PK.ms(1:length(xmc),i)=ms*PK.PixelSize;
            PK.Ms(1:length(xmc),i)=Ms*PK.PixelSize;
            PK.xmc(1:length(xmc),i)=xmc/PK.framerate;
            PK.xMc(1:length(xMc),i)=xMc/PK.framerate;
            PK.xmr(1:length(xmr),i)=xmr/PK.framerate;
            PK.xMr(1:length(xMr),i)=xMr/PK.framerate;
%             PK.locc(1:length(locc),i)=locc/PK.SamplingFrequency;
%             PK.locr(1:length(locr),i)=locr/PK.SamplingFrequency;
            
           

            
            
        end
        
        
        
     
        
        
        
    
        function PK=Save(PK,varargin)
            results_pathname=varargin{1};
            
        Tf= array2table(zeros(PK.number_cells,38));
        Tf.Properties.VariableNames = {'Period','Period_std','ascending_time','ascending_time_std','decay_time','decay_time_std',...
    'decay_time_90','decay_time_90_std','decay_time_70','decay_time_70_std','decay_time_50','decay_time_50_std','decay_time_30','decay_time_30_std','decay_time_20','decay_time_20_std',...
   'Tau_systole','Tau_systole_std','Baz_Tau_syst','Baz_Tau_syst_std','Tau_diast','Tau_diast_std','Baz_Tau_diast','Baz_Tau_diast_std','AUC','AUC_std','Pente_contraction','Pente_contraction_std','Pente_relax','Pente_relax_std',...
     'abs_amp_contraction','abs_amp_contraction_std','abs_amp_relax','abs_amp_relax_std','amp_max','amp_max_std','min','min_std'};

            MedT=nanmedian(diff(PK.posM,1,1),1);
            TabMedT=repmat(MedT,PK.ltab,1);
            Tf.Period=MedT';
            Tf.Period_std=nanstd(diff(PK.posM,1,1),1)';           
            Tf.ascending_time= nanmedian(PK.Tauc,1)';
            Tf.ascending_time_std= nanstd(PK.Tauc,1)';
            Tf.decay_time = nanmedian(PK.Taur,1)';
            Tf.decay_time_std= nanstd(PK.Taur,1)';
            Tf.Tau_systole = nanmedian(PK.Taus,1)';
            Tf.Tau_systole_std= nanstd(PK.Taus,1)';
            Tf.Baz_Tau_syst = nanmedian(PK.Taus./sqrt(TabMedT),1)';
            Tf.Baz_Tau_syst_std = nanstd(PK.Taus./sqrt(TabMedT),1)';
            Tf.Tau_diast = nanmedian(PK.Taud,1)';
            Tf.Tau_diast_std = nanstd(PK.Taud,1)';
            Tf.Baz_Tau_diast = nanmedian(PK.Taud./sqrt(TabMedT),1)';
            Tf.Baz_Tau_diast_std = nanstd(PK.Taud./sqrt(TabMedT),1)';
            Tf.AUC= nanmedian(PK.Aire,1)';
            Tf.AUC_std = nanstd(PK.Aire,1)';
            Tf.Pente_contraction = nanmedian(PK.pc,1)';
            Tf.Pente_contraction_std = nanstd(PK.pc,1)';
            Tf.Pente_relax= nanmedian(PK.pr,1)';
            Tf.Pente_relax_std = nanstd(PK.pr,1)';
            Tf.abs_amp_relax = nanmedian(PK.hr,1)';
            Tf.abs_amp_relax_std = nanstd(PK.hr,1)';
            Tf.abs_amp_contraction = nanmedian(PK.hc,1)';
            Tf.abs_amp_contraction_std = nanstd(PK.hc,1)';
            Tf.amp_max = nanmedian(PK.M,1)';
            Tf.amp_max_std = nanstd(PK.M,1)';
            Tf.min = nanmedian(PK.mmvg,1)';
            Tf.min_std = nanstd(PK.mmvg,1)';
            Tf.decay_time_20=nanmedian(PK.posper(:,1)-PK.posM,1)';
            Tf.decay_time_20_std=nanstd(PK.posper(:,1)-PK.posM,1)';
            Tf.decay_time_30=nanmedian(PK.posper(:,2)-PK.posM,1)';
            Tf.decay_time_30_std=nanstd(PK.posper(:,2)-PK.posM,1)';
            Tf.decay_time_50=nanmedian(PK.posper(:,3)-PK.posM,1)';
            Tf.decay_time_50_std=nanstd(PK.posper(:,3)-PK.posM)';
            Tf.decay_time_70=nanmedian(PK.posper(:,4)-PK.posM,1)';
            Tf.decay_time_70_std=nanstd(PK.posper(:,4)-PK.posM,1)';
            Tf.decay_time_90=nanmedian(PK.posper(:,5)-PK.posM,1)';
            Tf.decay_time_90_std=nanstd(PK.posper(:,5)-PK.posM,1)';
            
            
writetable(Tf,results_pathname,'WriteRowNames',true)           
            
        end
        
    end
    
end