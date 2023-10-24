classdef AnalysisPeaks < handle
    properties
        
        number_acquisitions
        number_cells
        vector_time
        matrix_rough_fluorescences
        vector_filtering_polynomial_order
        vector_filtering_frame_length
        matrix_filtered_fluorescences
        matrix_filtered_fluorescences_ori
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
        posMr
        ms
        Ms
        posper
        Mper
        ltab
        dd2
        smooth_signal
        base
        type
        nindpk
        basefit
        N
        f_smpks
        f_medpks
        f_multipks
        ind_smpks
        ind_medpks
        ind_multipks
        th_smpks
        th_medpks
        xmcr
        posmr
        mder
        Mder
        th_multi
        cut_freq
        type_bl
        bb
        bool_baselineref
        iter
        th
        Amea
        M2mea
        Tmea
        FPD
        xmcmea
        posMmea
        posmmea
        posM2mea
        mcmea
        Mmea
        mmea
        
    end
    methods
        function PK = AnalysisPeaks(varargin)
            p = inputParser;
            addRequired(p, 'Signal');
            addParameter(p, 'pkhgt',50)
            addOptional(p, 'PixelSize', 1);
            addOptional(p, 'SamplingFrequency', 1);
            addOptional(p, 'param_filter', 111);
            addOptional(p,'Pol_order',2)
            addOptional(p, 'Smoothness', 200);
            addOptional(p, 'proportion', 0.1);
            addOptional(p,'type',1)
            addOptional(p,'cut_freq',20)
            addOptional(p,'th_smpks',0.2)
            addOptional(p,'th_medpks',0.5)
            addOptional(p,'th_multi',1.5)
            addOptional(p,'baselinefit',0)
            addOptional(p,'bool_baselineref',0)
            parse(p, varargin{:});
            PK.iter=0;
            pol_length=p.Results.param_filter;
            sm=p.Results.Smoothness;
            th_smpks=p.Results.th_smpks;
            th_medpks=p.Results.th_medpks;
            th_multi=p.Results.th_multi;
            PK.bool_baselineref=p.Results.bool_baselineref;
            PK.bb=p.Results.baselinefit;
            PK.ltab=50000;
            prop=p.Results.proportion;
            PK.type=p.Results.type;
            pol_order=p.Results.Pol_order;
            PK.number_cells=size(p.Results.Signal,2)-1;
            time=p.Results.Signal(:,1);
            PK.PixelSize = p.Results.PixelSize;
            PK.cut_freq=p.Results.param_filter;
            
            %             PK.SamplingFrequency=p.Results.SamplingFrequency;
            %             PK.framerate = (PK.SamplingFrequency)/1000; %fps in ms
            
            ind=~isnan(time);
            
            PK.number_acquisitions=sum(ind);
            PK.vector_time=time(ind);
            PK.framerate=1/nanmedian(diff(PK.vector_time));
            PK.matrix_rough_fluorescences=p.Results.Signal(ind,2:(PK.number_cells+1));
            
            PK.vector_filtering_polynomial_order=pol_order*ones(1,PK.number_cells);
            PK.vector_filtering_frame_length=pol_length*ones(1,PK.number_cells);
            PK.matrix_filtered_fluorescences=zeros(PK.number_acquisitions,PK.number_cells);
            PK.smooth_signal=zeros(PK.number_acquisitions,PK.number_cells);
            PK.base=nan*ones(PK.number_acquisitions,PK.number_cells);
            PK.dd2=zeros(PK.number_acquisitions,PK.number_cells);
            PK.nindpk=true(round(PK.number_acquisitions),round(PK.number_cells));
            
            
            PK.sm=sm*ones(1,PK.number_cells);
            PK.th_smpks=th_smpks*ones(1,PK.number_cells);
            PK.th_medpks=th_medpks*ones(1,PK.number_cells);
            PK.th_multi=th_multi*ones(1,PK.number_cells);
            PK.prop=prop*ones(1,PK.number_cells);
            PK.xmr=nan*ones(PK.ltab,PK.number_cells);
            PK.ind_smpks=false(PK.ltab,PK.number_cells);
            PK.ind_medpks=false(PK.ltab,PK.number_cells);
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
        
        
        function PK=remove_base(PK,varargin)
            p=inputParser;
            addRequired(p,'i');
            parse(p,varargin{1:end});
            i=p.Results.i;
            
            
            PK.base(:,i)=nan*ones(PK.number_acquisitions,1);
            PK.base(PK.nindpk(:,i),i)= PK.matrix_filtered_fluorescences_ori(PK.nindpk(:,i),i);
            ind=~isnan(PK.base(:,i));
            switch PK.type_bl
                case 1
                    
                    ct=nanmean(PK.matrix_filtered_fluorescences_ori(ind,i));
                    PK.basefit(:,i)=ct*ones(length(PK.vector_time),1);
                case 2
                    titi=robustfit(PK.vector_time(ind),PK.matrix_filtered_fluorescences_ori(ind,i));
                    PK.basefit(:,i)=titi(2).*PK.vector_time+titi(1);
                case 3
                    ff=fit(PK.vector_time(ind),PK.matrix_filtered_fluorescences_ori(ind,i),'poly2');
                    PK.basefit(:,i)=ff(PK.vector_time);
                    
                case 4
                    ff=fit(PK.vector_time(ind),PK.matrix_filtered_fluorescences_ori(ind,i),'poly3');
                    PK.basefit(:,i)=ff(PK.vector_time);
            end
            
            
            
            
            
            PK.matrix_filtered_fluorescences(:,i)=PK.matrix_filtered_fluorescences_ori(:,i)-PK.basefit(:,i);
            PK.CalculateParameters(i);
            
        end
        
        function PK=Filter(PK,varargin)
            i=varargin{1};
            if PK.type<3
                PK.matrix_filtered_fluorescences(:,i)=sgolayfilt(PK.matrix_rough_fluorescences(:,i),PK.vector_filtering_polynomial_order(i),PK.vector_filtering_frame_length(i)); % in AU
                PK.matrix_filtered_fluorescences_ori(:,i)=PK.matrix_filtered_fluorescences(:,i);
                %  PK.matrix_filtered_fluorescences(:,i)=PK.matrix_rough_fluorescences(:,i);
            else
                F=fft(PK.matrix_rough_fluorescences(:,i));
                F(PK.cut_freq+1:end-PK.cut_freq)=0;
                PK.matrix_filtered_fluorescences(:,i)=real(ifft(F));
                PK.matrix_filtered_fluorescences_ori(:,i)=PK.matrix_filtered_fluorescences(:,i);
                
            end
            
        end
        
        
        
        
        function PK = CalculateParameters(PK,varargin)
            p=inputParser;
            addRequired(p,'i')
            parse(p, varargin{1:end});
            i=p.Results.i;
            dper=[0.8 0.7 0.5 0.3 0.05];
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
            PK.nindpk(:,i)=true(round(PK.number_acquisitions),1);
            PK.base(:,i)=nan*ones(round(PK.number_acquisitions),1);
            
            
            
            
            win=3;
            sl=PK.sm(i);
            %              Signal=PK.matrix_filtered_fluorescences(~isnan(PK.matrix_filtered_fluorescences(:,i)),i);
            Signal=PK.matrix_filtered_fluorescences(:,i);
            
            
            ss=smooth(Signal,sl,'loess');
            %                         ss=Signal;
%             figure;
%             hold on
%             plot(Signal,'linewidth',2)
            %                         plot(PK.vector_time,Signal)
            %                         plot(PK.matrix_rough_fluorescences(:,i),'color',[0.9 0.9 0.9])
            
            dd2=gradient(ss);
            PK.dd2(1:length(Signal),i)=dd2;
            PK.smooth_signal(1:length(Signal),i)=ss;
            dd=gradient(Signal);
            
            %                 if(sl~=80||prop~=0.05)
            %                                             figure
            %                                             hold on
            %    plot(dd)
            %                                             plot(dd2)
            %                                             pause
            %
            %                 end
            mm=max(dd2);
            nn=min(dd2);
            if PK.type<3
                PK.mder(i)=PK.prop(i)*nn;
                PK.Mder(i)=PK.prop(i)*mm;
            else
                PK.mder(i)=0.05*nn;
                PK.Mder(i)=0.05*mm;
            end
            
            [~, locct] = findpeaks(dd2,'MinPeakHeight' ,PK.Mder(i));
            [~, locrt] = findpeaks(-dd2,'MinPeakHeight' ,-PK.mder(i));
            
            clear maxc locc locr minr
            for ii=1:length(locct)
                [maxc(ii),locc(ii)]=max(dd(max(1,locct(ii)-win):min(locct(ii)+win,length(dd))));
                locc(ii)=max(1,locct(ii)-win)+locc(ii)-1;
            end
            for ii=1:length(locrt)
                [minr(ii),locr(ii)]=min(dd(max(1,locrt(ii)-win):min(locrt(ii)+win,length(dd))));
                locr(ii)=locr(ii)+max(1,locrt(ii)-win)-1;
            end
            %
            %                                                                              figure
            %                                                             hold on
            %                                                             plot(dd2)
            %                                                             plot(locc,maxc,'+r')
            %                                                             plot(locr,minr,'+g')
            %                                                             plot(sl,dd2(sl),'+y')
            %                                                             pause
            %
            
            [locr,ia,~]=unique(locr);
            minr=minr(ia);
            [locc,ia,~]=unique(locc);
            maxc=maxc(ia);
            
            
            %% for keeping the first peak in good cases
            
            
            [counts,ed]=histcounts(abs(dd2));
            T=otsuthresh(counts(2:end));
            
            %             figure
            %             plot(counts,'+')
            valt=T*max(ed);
            mini=std(abs(dd2(abs(dd2)<valt)));
            
            sl2=1;
            %             if abs(dd2(sl))<2*mini
            %                 locr=[sl2 locr];
            %                 minr=[dd2(sl2) minr];
            %             end
            %             windf=max(10,round(0.05*length(Signal)))
            windf=10;
            if abs(nanmean(dd2(sl:sl+windf)))<-PK.mder(i)/2
                locr=[sl2 locr];
                minr=[dd2(sl2) minr];
            end
            
            
            
            %% for keeping last peak in good cases
            
            %             if abs(dd2(end-sl))<2*mini
            %                 locc=[locc length(dd2)-sl2];
            %                 maxc=[maxc dd2(end-sl2)];
            %             end
            
            
            if abs(nanmean(dd2(end-sl-windf:end-sl)))<-PK.mder(i)/2||PK.type==3
                locc=[locc length(dd2)-sl2];
                maxc=[maxc dd2(end-sl2)];
            end
            %
            %
            
            %% calculate parameters for complete peaks
            
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
                            [M(count),posM(count)]=max(Signal(locc(ii):locr(pl)));
                            posM(count)=posM(count)+locc(ii)-1;
                            
                            
                            %                             [mmvg(count),vvg(ii)]=min(ss(locr(pl-1):locc(ii)));
                            [mmvg(count),vvg(ii)]=min(Signal(locr(pl-1):locc(ii)));
                            vvg(ii)=vvg(ii)+locr(pl-1)-1;
                            %                             [mmvd(count),vvd(ii)]=min(ss(locr(pl):locc(ii+1)));
                            [mmvd(count),vvd(ii)]=min(Signal(locr(pl):locc(ii+1)));
                            vvd(ii)=vvd(ii)+locr(pl)-1;
                            %                             win2=1;
                            %                             llgg=length(ss);
                            %
                            %                             mmvd(count)=nanmean(ss(max(1,vvd(ii)+win2):min(vvd(ii)+win2,llgg)));
                            %                             mmvg(count)=nanmean(ss(max(1,vvg(ii)-win2):min(vvg(ii)+win2,llgg)));
                            %                             plot(vvg(ii),mmvg(count),'+b')
                            %                             plot(vvd(ii),mmvd(count),'+r')
                            %                             plot(locr(pl-1),ss(locr(pl-1)),'+g')
                            %                             plot(locc(ii),ss(locc(ii)),'+m')
                            
                            
                            if PK.bool_baselineref&&PK.iter~=0
                                
                                mmvg(count)=0;
                                mmvd(count)=0;
                                
                            end
                            b=Signal(locc(ii))-maxc(ii)*locc(ii);
                            
                            xmc(count)=round((mmvg(count)-b)/maxc(ii));
                            xMc(count)=(M(count)-b)/maxc(ii);
                            b=Signal(locr(pl))-minr(pl)*locr(pl);
                            xMr(count)=(M(count)-b)/minr(pl);
                            xmr(count)=round((mmvd(count)-b)/minr(pl));
                            
                            %% correct the value of mmvd first calculated so that it fits on the curve
                            %                             mmvd(count)=Signal(round(xmr(count)));
                            %                                 %
                            %                             plot(xmc(count),mmvg(count),'+b')
                            %                             plot(xMc(count),M(count),'+y')
                            %                             plot(xmr(count),mmvd(count),'+g')
                            %                             plot(xMr(count),M(count),'+r')
                            
                            
                            
                            PK.pc(count,i)=maxc(ii)*PK.PixelSize*PK.framerate;
                            PK.hr(count,i)=(M(count)-mmvd(count))*PK.PixelSize;
                            PK.hc(count,i)=(M(count)-mmvg(count))*PK.PixelSize;
                            PK.M(count,i)=M(count)*PK.PixelSize;
                            PK.mmvg(count,i)=mmvg(count)*PK.PixelSize;
                            PK.mmvd(count,i)=mmvd(count)*PK.PixelSize;
                            
                            %                             [Ms(count),posM(count)]=max(Signal(round(max(1,xmc(count))):round(min(length(Signal),xmr(count)))));
                            %                             posM(count)=posM(count)+round(max(1,xmc(count)))-1;
                            
                            %                             plot(posM(count),Ms(count),'k+')
                            
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
            %             [ms(length(xMr)),posm(length(xMr))]=min(Signal(round(xMr(length(xMr))):end-sl));
            posm(length(xMr))=posm(length(xMr))+xMr(length(xMr))-1;
            %           plot(posm,ms,'+')
            %           plot(posM,Ms,'+')
            
            if PK.bool_baselineref&& PK.iter~=0
                ms=zeros(1,length(ms));
            end
            posper=round(repmat(xmr',1,5));
            for kk=1:length(xmr(:))
                for uu=1:5
                    % 5 is for 90% and 1 for 10%
                    dd=dper(uu)*(M(kk)-ms(kk));
                    tutu=find(diff(Signal(round(posM(kk)):round(posm(kk)))-ms(kk)>dd)==-1);
                    if ~isempty(tutu)
                        posper(kk,uu)=tutu(1);
                        posper(kk,uu)=posper(kk,uu)+round(posM(kk))-1;
                        Mper(kk,uu)=ms(kk)+dd;
                        
                        %                         Mper(kk,uu)=Signal(posper(kk,uu));
                        
                    end
                end
                
                pr(kk)=mean(dd2(max(posM(kk)+1,1):min(length(dd2),posper(kk,5)+1)));
                b=M(kk)-pr(kk)*posM(kk);
                xmrn(kk)=(mmvd(kk)-b)/pr(kk);
                
                if kk<length(xMr(:))
                    if xmrn(kk)<xmr(kk)||xmrn(kk)>xmc(kk+1)
                        xmrn(kk)=posm(kk);
                    end
                else
                    if xmrn(kk)<xmr(kk)||xmrn(kk)>length(Signal)
                        xmrn(kk)=posm(kk);
                    end
                end
                PK.Aire(kk,i)=sum(Signal(max(1,round(xmc(kk))):min(round(xmr(kk)),length(Signal))))/PK.framerate;
                xmrn=round(xmrn);
                PK.Tauc(kk,i)=PK.vector_time(posM(kk))-PK.vector_time(xmc(kk));
                PK.Taur(kk,i)=PK.vector_time(xmrn(kk))-PK.vector_time(posM(kk));
                PK.Taud(kk,i)=PK.vector_time(xmrn(kk))-PK.vector_time(xmc(kk));
                
                
                
                
                PK.nindpk(max(1,round(xmc(kk))):min(round(posper(kk,5)),length(Signal)),i)=0;
            end
            
            
            
            
            PK.nindpk(round(posm(end)):end,i)=0;
            PK.pr(1:length(xMr),i)=pr*PK.PixelSize*PK.framerate;
            
            PK.ms(1:length(xmc),i)=ms*PK.PixelSize;
            %            PK.Ms(1:length(xmc),i)=Ms*PK.PixelSize;
            PK.posper(1:size(posper,1),:,i)=reshape(PK.vector_time(posper(:)+1),size(posper));
            PK.Mper(1:size(posper,1),:,i)=reshape(Signal(posper(:)+1),size(posper));
            
            %             PK.posm(1:length(xmc),i)=(posm-1)/PK.framerate;
            %             PK.posmr(1:length(xmc),i)=posm;
            %             PK.posM(1:length(xmc),i)=(posM-1)/PK.framerate;
            %             PK.ms(1:length(xmc),i)=ms*PK.PixelSize;
            %             PK.Ms(1:length(xmc),i)=Ms*PK.PixelSize;
            %             PK.xmc(1:length(xmc),i)=(xmc-1)/PK.framerate;
            %             PK.xmcr(1:length(xmc),i)=xmc;
            %             PK.xMc(1:length(xMc),i)=(xMc-1)/PK.framerate;
            %             PK.xmr(1:length(xmr),i)=(xmr-1)/PK.framerate;
            %             PK.xMr(1:length(xMr),i)=(xMr-1)/PK.framerate;
            
            
            PK.posm(1:length(xmc),i)=PK.vector_time(round(posm));
            PK.posmr(1:length(xmc),i)=round(posm);
            PK.posM(1:length(xmc),i)=PK.vector_time(posM);
            PK.posMr(1:length(xmc),i)=round(posM);
            PK.xmc(1:length(xmc),i)=PK.vector_time(round(xmc));
            PK.xmcr(1:length(xmc),i)=round(xmc);
            %         PK.xMc(1:length(xMc),i)=PK.vector_time(round(xMc));
            PK.xmr(1:length(xmr),i)=PK.vector_time(round(xmrn));
            %      PK.xMr(1:length(xMr),i)=PK.vector_time(round(xMr));
            
            
            %             PK.locc(1:length(locc),i)=locc/PK.SamplingFrequency;
            %             PK.locr(1:length(locr),i)=locr/PK.SamplingFrequency;
            
            %             plot(PK.posM(:,i),PK.Ms(:,i),'+k')
            PK.N(i)=length(PK.posM(~isnan(PK.posM(:,i)),i));
            PK.iter=PK.iter+1;
            
        end
        
        function PK=statpks(PK,varargin)
            
            p=inputParser;
            addRequired(p,'i')
            parse(p, varargin{1:end});
            i=p.Results.i;
            
            
            PK.ind_smpks(:,i)=false(PK.ltab,1);
            PK.ind_medpks(:,i)=false(PK.ltab,1);
            PK.ind_multipks(:,i)=false(PK.ltab,1);
            contsmpks=0;
            contmultipks=0;
            contmedpks=0;
            MM=nanmax(cat(1,PK.hr(:,i),PK.hc(:,i)));
            titi=cat(1,PK.mmvg(~isnan(PK.mmvg(:,i)),i),PK.mmvd(~isnan(PK.mmvd(:,i)),i));
            bb=quantile(titi,0.4);
            %             bb=nanmedian(cat(1,PK.mmvg(:,i),PK.mmvd(:,i)));
            tutu=cat(1,PK.hr(~isnan(PK.hr(:,i)),i),PK.hc(~isnan(PK.hc(:,i)),i));
            Q=quantile(tutu,0.9);
            Medpk=PK.th_multi(i)*Q;
            %             bb_std=PK.th_multi(i)*nanstd(cat(1,PK.mmvg(:,i),PK.mmvd(:,i)));
            
            for ii=1:length(PK.posM(~isnan(PK.posM(:,i)),i))
                
                if max(PK.hr(ii,i),PK.hc(ii,i))<MM*PK.th_smpks(i)
                    contsmpks=contsmpks+1;
                    PK.ind_smpks(ii,i)=1;
                end
                if max(PK.hr(ii,i),PK.hc(ii,i))<MM*PK.th_medpks(i)&& max(PK.hr(ii,i),PK.hc(ii,i))>MM*PK.th_smpks(i)
                    contmedpks=contmedpks+1;
                    PK.ind_medpks(ii,i)=1;
                end
                
                
                if((PK.mmvg(ii,i))>(bb+Medpk))||((PK.mmvd(ii,i))>(bb+Medpk))
                    contmultipks=contmultipks+1;
                    PK.ind_multipks(ii,i)=1;
                end
                
            end
            
            
            PK.f_smpks(i)=contsmpks/PK.N(i)*100;
            PK.f_medpks(i)=contmedpks/PK.N(i)*100;
            PK.f_multipks(i)=contmultipks/PK.N(i)*100;
            
        end
        
        function PK=Analyse_MEA(PK,varargin)
            p=inputParser;
            addRequired(p,'i')
            parse(p, varargin{1:end});
            i=p.Results.i;
            Signal=PK.matrix_filtered_fluorescences(:,i);
            mm=max(Signal);
            
            PK.posmmea(:,i)=nan*ones(PK.ltab,1);
            PK.posMmea(:,i)=nan*ones(PK.ltab,1);
            PK.mmea(:,i)=nan*ones(PK.ltab,1);
            PK.Mmea(:,i)=nan*ones(PK.ltab,1);
            PK.xmcmea(:,i)=nan*ones(PK.ltab,1);
            PK.mcmea(:,i)=nan*ones(PK.ltab,1);
            PK.posM2mea(:,i)=nan*ones(PK.ltab,1);
            PK.M2mea(:,i)=nan*ones(PK.ltab,1);
            PK.Amea(:,i)=nan*ones(PK.ltab,1);
            PK.FPD(:,i)=nan*ones(PK.ltab,1);
            PK.Tmea(:,i)=nan*ones(PK.ltab,1);

            
            
            
            PK.th(i)=PK.prop(i)*mm;
            
            [M,posM] = findpeaks(Signal,'MinPeakHeight' ,PK.th(i));
      
            %             figure
            %             hold on
            %             plot(Signal)
            
            %% find parameters for complete peaks
            count=1;
            %             try
            for uu=1: length(posM)-1
           
                    ind=posM(uu):posM(uu+1);
    
                [m(uu),posm(uu)]=min(Signal(ind));
                posm(uu)=posm(uu)+posM(uu)-1;
                if uu<length(posM)
             
%                     aa=PK.posMr==posM(uu+1);
%                     delta=(posM(uu+1)-PK.xmcr(aa))/5;
%                     xmcmea(uu)=PK.xmcr(aa);
%                     ind2=posm(uu):PK.xmcr(aa)-delta;
%             
%                     [M2(uu),posM2(uu)]=max(Signal(ind2));
%                     posM2(uu)=posM2(uu)+posm(uu)-1;
            % other simpler method
                ind2=posm(uu):posM(uu+1);
                [Mt,post]=findpeaks(Signal(ind2),'SortStr','descend');
                M2(uu)=Mt(1);
                posM2(uu)=post(1)+posm(uu)-1;


                end
                
                
            end
            
            %                plot(posM,M,'+r')
            %             plot(posm,m,'+y')
            %             plot(posM2,M2,'+g')
            %             catch
            %             end
  

            PK.Amea(1:length(M)-1,i)=M(1:end-1)'-m;
            PK.M2mea(1:length(M2),i)=M2;
            PK.Tmea(1:length(posM)-1,i)=(diff(posM))/PK.framerate;
            PK.FPD(1:length(posm),i)=(posM2-posm)/PK.framerate;
%             PK.xmcmea(1:length(xmcmea),i)=PK.vector_time(xmcmea);
%             PK.mcmea(1:length(xmcmea),i)=Signal(xmcmea);
            PK.posMmea(1:length(posM),i)=PK.vector_time(posM);
            PK.posmmea(1:length(posm),i)=PK.vector_time(posm);
            PK.posM2mea(1:length(posM2),i)=PK.vector_time(posM2);
            PK.Mmea(1:length(M),i)=M;
            PK.mmea(1:length(m),i)=m;
            
        end
        
        
        
        
        
        function PK=Save(PK,varargin)
            results_pathname=varargin{1};
            
            switch PK.type
                case 2
                    Tf= array2table(zeros(PK.number_cells,35));
                    
                    Tf.Properties.VariableNames = {'N_pks','Period','Period_std','ascending_time','ascending_time_std','decay_time','decay_time_std',...
                        'decay_time_90','decay_time_90_std','decay_time_70','decay_time_70_std','decay_time_50','decay_time_50_std','decay_time_30','decay_time_30_std','decay_time_20','decay_time_20_std',...
                        'Tau_d','Tau_d_std','Baz_Tau_d','Baz_Tau_d_std','AUC','AUC_std','Pente_contraction','Pente_contraction_std','Pente_relax','Pente_relax_std',...
                        'abs_amp_contraction','abs_amp_contraction_std','abs_amp_relax','abs_amp_relax_std','amp_max','amp_max_std','min','min_std'};
                    
                case 1
                    
                    
                    Tf= array2table(zeros(PK.number_cells,38));
                    Tf.Properties.VariableNames = {'Period','Period_std','ascending_time','ascending_time_std','decay_time','decay_time_std',...
                        'decay_time_90','decay_time_90_std','decay_time_70','decay_time_70_std','decay_time_50','decay_time_50_std','decay_time_30','decay_time_30_std','decay_time_20','decay_time_20_std',...
                        'Tau_d','Tau_d_std','Baz_Tau_d','Baz_Tau_d_std','AUC','AUC_std','Pente_contraction','Pente_contraction_std','Pente_relax','Pente_relax_std',...
                        'abs_amp_contraction','abs_amp_contraction_std','abs_amp_relax','abs_amp_relax_std','amp_max','amp_max_std','min','min_std',...
                        'N_pks','f_smpks','f_medpks','f_multipks'};
                    
                case 3
                    Tf= array2table(zeros(PK.number_cells,6));
                    Tf.Properties.VariableNames = {'Period','Period_std','FP_Duration','FP_Duration_std','Amplitude','Amplitude_std'};
                    
                    
            end
            MedT=NaN*ones(1,PK.number_cells);
            StdT=NaN*ones(1,PK.number_cells);
            for uu=1:PK.number_cells
                MedT(uu)=median(diff(PK.posM(~isnan(PK.posM(:,uu)),uu)));
                StdT(uu)=std(diff(PK.posM(~isnan(PK.posM(:,uu)),uu)));
            end
            TabMedT=repmat(MedT,PK.ltab,1);
            
            
            if PK.type<3
                Tf.Period=MedT';
                Tf.Period_std=StdT';
                Tf.ascending_time= nanmean(PK.Tauc,1)';
                Tf.ascending_time_std= nanstd(PK.Tauc,1)';
                Tf.decay_time = nanmean(PK.Taur,1)';
                Tf.decay_time_std= nanstd(PK.Taur,1)';
                %             Tf.Tau_systole = nanmean(PK.Taus,1)';
                %             Tf.Tau_systole_std= nanstd(PK.Taus,1)';
                %             Tf.Baz_Tau_syst = nanmean(PK.Taus./sqrt(TabMedT),1)';
                %             Tf.Baz_Tau_syst_std = nanstd(PK.Taus./sqrt(TabMedT),1)';
                Tf.Tau_d = nanmean(PK.Taud,1)';
                Tf.Tau_d_std = nanstd(PK.Taud,1)';
                Tf.Baz_Tau_d = nanmean(PK.Taud./sqrt(TabMedT),1)';
                Tf.Baz_Tau_d_std = nanstd(PK.Taud./sqrt(TabMedT),1)';
                Tf.AUC= nanmean(PK.Aire,1)';
                Tf.AUC_std = nanstd(PK.Aire,1)';
                Tf.Pente_contraction = nanmean(PK.pc,1)';
                Tf.Pente_contraction_std = nanstd(PK.pc,1)';
                Tf.Pente_relax= nanmean(PK.pr,1)';
                Tf.Pente_relax_std = nanstd(PK.pr,1)';
                Tf.abs_amp_relax = nanmean(PK.hr,1)';
                Tf.abs_amp_relax_std = nanstd(PK.hr,1)';
                Tf.abs_amp_contraction = nanmean(PK.hc,1)';
                Tf.abs_amp_contraction_std = nanstd(PK.hc,1)';
                Tf.amp_max = nanmean(PK.M,1)';
                Tf.amp_max_std = nanstd(PK.M,1)';
                Tf.min = nanmean(PK.mmvg,1)';
                Tf.min_std = nanstd(PK.mmvg,1)';
                Tf.decay_time_20=nanmean(squeeze(PK.posper(:,1,:))-PK.posM,1)';
                Tf.decay_time_20_std=nanstd(squeeze(PK.posper(:,1,:))-PK.posM,1)';
                Tf.decay_time_30=nanmean(squeeze(PK.posper(:,2,:))-PK.posM,1)';
                Tf.decay_time_30_std=nanstd(squeeze(PK.posper(:,2,:))-PK.posM,1)';
                Tf.decay_time_50=nanmean(squeeze(PK.posper(:,3,:))-PK.posM,1)';
                Tf.decay_time_50_std=nanstd(squeeze(PK.posper(:,3,:))-PK.posM)';
                Tf.decay_time_70=nanmean(squeeze(PK.posper(:,4,:))-PK.posM,1)';
                Tf.decay_time_70_std=nanstd(squeeze(PK.posper(:,4,:))-PK.posM,1)';
                Tf.decay_time_90=nanmean(squeeze(PK.posper(:,5,:))-PK.posM,1)';
                Tf.decay_time_90_std=nanstd(squeeze(PK.posper(:,5,:))-PK.posM,1)';
                Tf.N_pks=PK.N';
                if(PK.type==1)
                    Tf.f_smpks=PK.f_smpks';
                    Tf.f_medpks=PK.f_medpks';
                    Tf.f_multipks=PK.f_multipks';
                end
            end
            if PK.type==3
                Tf.Period=nanmean(PK.Tmea,1)';
                Tf.Period_std=nanstd(PK.Tmea,1)';
                Tf.Amplitude=nanmean(PK.Amea,1)';
                Tf.Amplitude_std=nanstd(PK.Amea,1)';
                Tf.FP_Duration=nanmean(PK.FPD,1)';
                Tf.FP_Duration_std=nanstd(PK.FPD,1)';
                
                
            end
            
            
            Tfa = table2array(Tf);
            Tff = array2table(Tfa.');
            Tff.Properties.RowNames = Tf.Properties.VariableNames;
            writetable(Tff,results_pathname,'WriteRowNames',true)
            
        end
        
    end
    
end