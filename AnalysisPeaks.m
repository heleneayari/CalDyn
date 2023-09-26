classdef AnalysisPeaks < handle
    properties
        
        number_acquisitions
        number_cells
        vector_time
        matrix_rough_fluorescences
        vector_filtering_polynomial_order
        vector_filtering_frame_length
        matrix_filtered_fluorescences
        vector_threshold_peaks_detection
        vector_number_peaks
        vector_number_valleys
        matrix_time_peaks
        matrix_time_valleys
        matrix_fluorescences_peaks
        matrix_fluorescences_valleys
        vector_number_complete_peaks
        vector_number_complete_valleys
        matrix_time_complete_peaks
        matrix_time_complete_valleys
        matrix_fluorescences_complete_peaks
        matrix_fluorescences_complete_valleys
        matrix_left_durations_complete_peaks
        matrix_right_durations_complete_peaks
        matrix_left_amplitudes_complete_peaks
        matrix_left_normalized_complete_peaks
        matrix_right_amplitudes_complete_peaks
        matrix_right_normalized_amplitudes_complete_peaks
        matrix_left_velocities_complete_peaks
        matrix_right_velocities_complete_peaks
        vector_frequency_complete_peaks
        vector_frequency_complete_valleys
        matrix_effective_areas_complete_peaks
        vector_maximum_left_amplitude_complete_peaks
        vector_maximum_left_normalized_amplitude_complete_peaks
        vector_maximum_right_amplitude_complete_peaks
        vector_maximum_right_normalized_amplitude_complete_peaks
        vector_threshold_small_medium_peaks
        vector_threshold_medium_large_peaks
        vector_number_small_complete_peaks
        vector_frequency_small_complete_peaks
        matrix_numbers_small_complete_peaks
        matrix_left_durations_small_complete_peaks
        matrix_right_durations_small_complete_peaks
        matrix_left_amplitudes_small_complete_peaks
        matrix_left_normalized_amplitudes_small_complete_peaks
        matrix_right_amplitudes_small_complete_peaks
        matrix_right_normalized_amplitudes_small_complete_peaks
        matrix_left_velocities_small_complete_peaks
        matrix_right_velocities_small_complete_peaks
        matrix_areas_small_complete_peaks
        vector_number_medium_complete_peaks
        vector_frequency_medium_complete_peaks
        matrix_numbers_medium_complete_peaks
        matrix_left_durations_medium_complete_peaks
        matrix_right_durations_medium_complete_peaks
        matrix_left_amplitudes_medium_complete_peaks
        matrix_left_normalized_amplitudes_medium_complete_peaks
        matrix_right_amplitudes_medium_complete_peaks
        matrix_right_normalized_amplitudes_medium_complete_peaks
        matrix_left_velocities_medium_complete_peaks
        matrix_right_velocities_medium_complete_peaks
        matrix_areas_medium_complete_peaks
        vector_number_large_complete_peaks
        vector_frequency_large_complete_peaks
        matrix_numbers_large_complete_peaks
        matrix_left_durations_large_complete_peaks
        matrix_right_durations_large_complete_peaks
        matrix_left_amplitudes_large_complete_peaks
        matrix_left_normalized_amplitudes_large_complete_peaks
        matrix_right_amplitudes_large_complete_peaks
        matrix_right_normalized_amplitudes_large_complete_peaks
        matrix_left_velocities_large_complete_peaks
        matrix_right_velocities_large_complete_peaks
        matrix_areas_large_complete_peaks
        vector_fluorescences_peaks
        vector_time_peaks
        vector_fluorescences_valleys
        vector_time_valleys
        vector_time_complete_peaks
        vector_time_complete_valleys
        vector_fluorescences_complete_peaks
        vector_fluorescences_complete_valleys
        vector_left_durations_complete_peaks
        vector_right_durations_complete_peaks
        vector_left_amplitudes_complete_peaks
        vector_left_normalized_amplitudes_complete_peaks
        vector_right_amplitudes_complete_peaks
        vector_right_normalized_amplitudes_complete_peaks
        vector_left_velocities_complete_peaks
        vector_right_velocities_complete_peaks
        vector_numbers_small_complete_peaks
        vector_left_durations_small_complete_peaks
        vector_right_durations_small_complete_peaks
        vector_left_amplitudes_small_complete_peaks
        vector_left_normalized_amplitudes_small_complete_peaks
        vector_right_amplitudes_small_complete_peaks
        vector_right_normalized_amplitudes_small_complete_peaks
        vector_left_velocities_small_complete_peaks
        vector_right_velocities_small_complete_peaks
        vector_areas_small_complete_peaks
        vector_numbers_medium_complete_peaks
        vector_left_durations_medium_complete_peaks
        vector_right_durations_medium_complete_peaks
        vector_left_amplitudes_medium_complete_peaks
        vector_left_normalized_amplitudes_medium_complete_peaks
        vector_right_amplitudes_medium_complete_peaks
        vector_right_normalized_amplitudes_medium_complete_peaks
        vector_left_velocities_medium_complete_peaks
        vector_right_velocities_medium_complete_peaks
        vector_effective_areas_complete_peaks
        vector_areas_medium_complete_peaks
        vector_numbers_large_complete_peaks
        vector_left_durations_large_complete_peaks
        vector_right_durations_large_complete_peaks
        vector_left_amplitudes_large_complete_peaks
        vector_left_normalized_amplitudes_large_complete_peaks
        vector_right_amplitudes_large_complete_peaks
        vector_right_normalized_amplitudes_large_complete_peaks
        vector_left_velocities_large_complete_peaks
        vector_right_velocities_large_complete_peaks
        vector_total_areas_complete_peaks
        vector_areas_large_complete_peaks
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
        
        
    end
    methods
        function PK = AnalysisPeaks(varargin)
            p = inputParser;
            addRequired(p, 'Signal');
            addParameter(p, 'pkhgt',50)
            addOptional(p, 'PixelSize', 1);
            addOptional(p, 'SamplingFrequency', 1);
            parse(p, varargin{:});
            pkh=p.Results.pkhgt;
            PK.PixelSize = p.Results.PixelSize;
            PK.SamplingFrequency=p.Results.SamplingFrequency;
            PK.framerate = (PK.SamplingFrequency)/1000; %fps in ms
            PK.number_acquisitions=size(p.Results.Signal,1);
            PK.number_cells=size(p.Results.Signal,2)-1;
            PK.vector_time=p.Results.Signal(:,1);
            PK.matrix_rough_fluorescences=p.Results.Signal(:,2:(PK.number_cells+1));
            PK.vector_filtering_polynomial_order=2*ones(1,PK.number_cells);
            PK.vector_filtering_frame_length=11*ones(1,PK.number_cells);
            PK.matrix_filtered_fluorescences=zeros(PK.number_acquisitions,PK.number_cells);
            PK.vector_threshold_peaks_detection=pkh*ones(1,PK.number_cells);
            PK.vector_number_peaks=zeros(1,PK.number_cells);
            PK.vector_number_valleys=zeros(1,PK.number_cells);
            PK.matrix_time_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_time_valleys=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_fluorescences_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_fluorescences_valleys=zeros(PK.number_acquisitions,PK.number_cells);
            PK.vector_number_complete_peaks=zeros(1,PK.number_cells);
            PK.vector_number_complete_valleys=zeros(1,PK.number_cells);
            PK.matrix_time_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_time_complete_valleys=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_fluorescences_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_fluorescences_complete_valleys=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_left_durations_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_right_durations_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_left_amplitudes_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_left_normalized_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_right_amplitudes_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_right_normalized_amplitudes_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_left_velocities_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_right_velocities_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.vector_frequency_complete_peaks=zeros(1,PK.number_cells);
            PK.vector_frequency_complete_valleys=zeros(1,PK.number_cells);
            PK.matrix_effective_areas_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.vector_maximum_left_amplitude_complete_peaks=zeros(1,PK.number_cells);
            PK.vector_maximum_left_normalized_amplitude_complete_peaks=zeros(1,PK.number_cells);
            PK.vector_maximum_right_amplitude_complete_peaks=zeros(1,PK.number_cells);
            PK.vector_maximum_right_normalized_amplitude_complete_peaks=zeros(1,PK.number_cells);
            PK.vector_threshold_small_medium_peaks=0.2*ones(1,PK.number_cells);
            PK.vector_threshold_medium_large_peaks=0.5*ones(1,PK.number_cells);
            PK.vector_number_small_complete_peaks=zeros(1,PK.number_cells);
            PK.vector_frequency_small_complete_peaks=zeros(1,PK.number_cells);
            PK.matrix_numbers_small_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_left_durations_small_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_right_durations_small_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_left_amplitudes_small_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_left_normalized_amplitudes_small_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_right_amplitudes_small_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_right_normalized_amplitudes_small_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_left_velocities_small_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_right_velocities_small_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_areas_small_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.vector_number_medium_complete_peaks=zeros(1,PK.number_cells);
            PK.vector_frequency_medium_complete_peaks=zeros(1,PK.number_cells);
            PK.matrix_numbers_medium_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_left_durations_medium_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_right_durations_medium_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_left_amplitudes_medium_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_left_normalized_amplitudes_medium_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_right_amplitudes_medium_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_right_normalized_amplitudes_medium_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_left_velocities_medium_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_right_velocities_medium_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_areas_medium_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.vector_number_large_complete_peaks=zeros(1,PK.number_cells); %%%%%%
            PK.vector_frequency_large_complete_peaks=zeros(1,PK.number_cells); %%%%%%
            PK.matrix_numbers_large_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_left_durations_large_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_right_durations_large_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_left_amplitudes_large_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_left_normalized_amplitudes_large_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_right_amplitudes_large_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_right_normalized_amplitudes_large_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_left_velocities_large_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_right_velocities_large_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            PK.matrix_areas_large_complete_peaks=zeros(PK.number_acquisitions,PK.number_cells);
            
            PK.sm=20*ones(1,PK.number_cells);
            PK.prop=0.2*ones(1,PK.number_cells);
            PK.xmr=nan*ones(100,PK.number_cells);
            PK.xMr=nan*ones(100,PK.number_cells);
            PK.xmc=nan*ones(100,PK.number_cells);
            PK.xMc=nan*ones(100,PK.number_cells);
            PK.mmvg=nan*ones(100,PK.number_cells);
            PK.mmvd=nan*ones(100,PK.number_cells);
            PK.M=nan*ones(100,PK.number_cells);
            PK.Tauc=nan*ones(100,PK.number_cells);
            PK.Taur=nan*ones(100,PK.number_cells);
            PK.Taud=nan*ones(100,PK.number_cells);
            PK.Taus=nan*ones(100,PK.number_cells);
            PK.Tauc=nan*ones(100,PK.number_cells);
            PK.Tauc=nan*ones(100,PK.number_cells);
            PK.pc=nan*ones(100,PK.number_cells);
            PK.pr=nan*ones(100,PK.number_cells);
            PK.hc=nan*ones(100,PK.number_cells);
            PK.hr=nan*ones(100,PK.number_cells);
            PK.Aire=nan*ones(100,PK.number_cells);
            PK.locc=nan*ones(100,PK.number_cells);
            PK.locr=nan*ones(100,PK.number_cells);
            PK.posm=nan*ones(100,PK.number_cells);
            PK.posM=nan*ones(100,PK.number_cells);
            PK.ms=nan*ones(100,PK.number_cells);
            PK.Ms=nan*ones(100,PK.number_cells);
            PK.posper=nan*ones(100,5,PK.number_cells);
            PK.Mper=nan*ones(100,5,PK.number_cells);     
            
            
        end
        
        
        function PK=Filter(PK,varargin)
            i=varargin{1};
            PK.matrix_filtered_fluorescences(:,i)=sgolayfilt(PK.matrix_rough_fluorescences(:,i),PK.vector_filtering_polynomial_order(i),PK.vector_filtering_frame_length(i)); % in AU
%             PK.matrix_filtered_fluorescences(:,i)=PK.matrix_rough_fluorescences(:,i);
        end
        
        function PK=Thresh(PK,varargin)
            i=varargin{1};
     
            [PK.vector_fluorescences_peaks,PK.vector_time_peaks]=findpeaks(PK.matrix_filtered_fluorescences(:,i),PK.vector_time,'MinPeakProminence',PK.vector_threshold_peaks_detection(i));
            [PK.vector_fluorescences_valleys,PK.vector_time_valleys]=findpeaks(-PK.matrix_filtered_fluorescences(:,i),PK.vector_time,'MinPeakProminence',PK.vector_threshold_peaks_detection(i));
            if length(PK.vector_time_peaks)==1 && isempty(PK.vector_time_valleys)
                PK.vector_fluorescences_valleys(1,1)=max(-PK.matrix_filtered_fluorescences(1:find(PK.vector_time==PK.vector_time_peaks),1));
                PK.vector_time_valleys(1,1)=PK.vector_time(-PK.matrix_filtered_fluorescences(:,i)==PK.vector_fluorescences_valleys(1,1));
                PK.vector_fluorescences_valleys(2,1)=max(-PK.matrix_filtered_fluorescences(find(PK.vector_time==PK.vector_time_peaks):length(PK.vector_time),1));
                PK.vector_time_valleys(2,1)=PK.vector_time(-PK.matrix_filtered_fluorescences(:,i)==PK.vector_fluorescences_valleys(2,1));
            end
            PK.vector_fluorescences_valleys=-PK.vector_fluorescences_valleys; % in AU
            PK.vector_number_peaks(i)=length(PK.vector_time_peaks);
            PK.vector_number_valleys(i)=length(PK.vector_time_valleys);
            PK.matrix_time_peaks(1:PK.vector_number_peaks(i),i)=PK.vector_time_peaks; % in s
            PK.matrix_time_valleys(1:PK.vector_number_valleys(i),i)=PK.vector_time_valleys; % in s
            PK.matrix_fluorescences_peaks(1:PK.vector_number_peaks(i),i)=PK.vector_fluorescences_peaks; % in AU
            PK.matrix_fluorescences_valleys(1:PK.vector_number_valleys(i),i)=PK.vector_fluorescences_valleys; % in AU
            
%             close(figure(2));
%             figure(2)
%             plot(PK.vector_time,PK.matrix_filtered_fluorescences(:,i),'-k');
%             grid on
%             title('Filtered signal','interpreter','latex')
%             xlabel('Time (s)','interpreter','latex');
%             ylabel('Fluorescence (AU)','interpreter','latex');
%             hold on
%             plot(PK.vector_time_peaks,PK.vector_fluorescences_peaks,'+r');
%             hold on
%             plot(PK.vector_time_valleys,PK.vector_fluorescences_valleys,'xb');
%             legend('Experimental data (filtered)','Peaks','Valleys');
            
        end
        function PK=Calculate(PK,varargin)
            i=varargin{1};

            
            
            
            
            
            
            
            
            
            
            
            
            
            
            if PK.vector_time_valleys(1)<PK.vector_time_peaks(1)
                if PK.vector_number_peaks(i)<PK.vector_number_valleys(i)
                    PK.vector_number_complete_peaks(i)=PK.vector_number_peaks(i);
                    PK.vector_number_complete_valleys(i)=PK.vector_number_valleys(i);
                    PK.vector_time_complete_peaks=PK.vector_time_peaks; % in s
                    PK.vector_time_complete_valleys=PK.vector_time_valleys; % in s
                    PK.vector_fluorescences_complete_peaks=PK.vector_fluorescences_peaks; % in AU
                    PK.vector_fluorescences_complete_valleys=PK.vector_fluorescences_valleys; % in AU
                else
                    PK.vector_number_complete_peaks(i)=PK.vector_number_peaks(i)-1;
                    PK.vector_number_complete_valleys(i)=PK.vector_number_valleys(i);
                    PK.vector_time_complete_peaks=PK.vector_time_peaks(1:PK.vector_number_complete_peaks(i)); % in s
                    PK.vector_time_complete_valleys=PK.vector_time_valleys; % in s
                    PK.vector_fluorescences_complete_peaks=PK.vector_fluorescences_peaks(1:PK.vector_number_complete_peaks(i)); % in AU
                    PK.vector_fluorescences_complete_valleys=PK.vector_fluorescences_valleys; % in AU
                end
            else
                if PK.vector_number_peaks(i)==PK.vector_number_valleys(i)
                    PK.vector_number_complete_peaks(i)=PK.vector_number_peaks(i)-1;
                    PK.vector_number_complete_valleys(i)=PK.vector_number_valleys(i);
                    PK.vector_time_complete_peaks=PK.vector_time_peaks(2:(PK.vector_number_complete_peaks(i)+1)); % in s
                    PK.vector_time_complete_valleys=PK.vector_time_valleys; % in s
                    PK.vector_fluorescences_complete_peaks=PK.vector_fluorescences_peaks(2:(PK.vector_number_complete_peaks(i)+1)); % in AU
                    PK.vector_fluorescences_complete_valleys=PK.vector_fluorescences_valleys; % in AU
                else
                    PK.vector_number_complete_peaks(i)=PK.vector_number_peaks(i)-2;
                    PK.vector_number_complete_valleys(i)=PK.vector_number_valleys(i);
                    PK.vector_time_complete_peaks=PK.vector_time_peaks(2:(PK.vector_number_complete_peaks(i)+1)); % in s
                    PK.vector_time_complete_valleys=PK.vector_time_valleys; % in s
                    PK.vector_fluorescences_complete_peaks=PK.vector_fluorescences_peaks(2:(PK.vector_number_complete_peaks(i)+1)); % in AU
                    PK.vector_fluorescences_complete_valleys=PK.vector_fluorescences_valleys; % in AU
                end
            end
            PK.vector_left_durations_complete_peaks=PK.vector_time_complete_peaks-PK.vector_time_complete_valleys(1:(PK.vector_number_complete_valleys(i)-1)); % in s
            PK.vector_right_durations_complete_peaks=-(PK.vector_time_complete_peaks-PK.vector_time_complete_valleys(2:PK.vector_number_complete_valleys(i))); % in s
            PK.vector_left_amplitudes_complete_peaks=PK.vector_fluorescences_complete_peaks-PK.vector_fluorescences_complete_valleys(1:(PK.vector_number_complete_valleys(i)-1)); % in AU
            PK.vector_left_normalized_amplitudes_complete_peaks=PK.vector_left_amplitudes_complete_peaks./PK.vector_fluorescences_complete_valleys(1:(PK.vector_number_complete_valleys(i)-1));
            PK.vector_right_amplitudes_complete_peaks=PK.vector_fluorescences_complete_peaks-PK.vector_fluorescences_complete_valleys(2:PK.vector_number_complete_valleys(i)); % in AU
            PK.vector_right_normalized_amplitudes_complete_peaks=PK.vector_right_amplitudes_complete_peaks./PK.vector_fluorescences_complete_valleys(2:PK.vector_number_complete_valleys(i));
            PK.vector_left_velocities_complete_peaks=PK.vector_left_amplitudes_complete_peaks./PK.vector_left_durations_complete_peaks; % in AU/s
            PK.vector_right_velocities_complete_peaks=PK.vector_right_amplitudes_complete_peaks./PK.vector_right_durations_complete_peaks; % in AU/s
            PK.matrix_time_complete_peaks(1:PK.vector_number_complete_peaks(i),i)=PK.vector_time_complete_peaks; % in s
            PK.matrix_time_complete_valleys(1:PK.vector_number_complete_valleys(i),i)=PK.vector_time_complete_valleys; % in s
            PK.matrix_fluorescences_complete_peaks(1:PK.vector_number_complete_peaks(i),i)=PK.vector_fluorescences_complete_peaks; % in AU
            PK.matrix_fluorescences_complete_valleys(1:PK.vector_number_complete_valleys(i),i)=PK.vector_fluorescences_complete_valleys; % in AU
            PK.matrix_left_durations_complete_peaks(1:PK.vector_number_complete_peaks(i),i)=PK.vector_left_durations_complete_peaks; % in s
            PK.matrix_right_durations_complete_peaks(1:PK.vector_number_complete_peaks(i),i)=PK.vector_right_durations_complete_peaks; % in s
            PK.matrix_left_amplitudes_complete_peaks(1:PK.vector_number_complete_peaks(i),i)=PK.vector_left_amplitudes_complete_peaks; % in AU
            PK.matrix_left_normalized_complete_peaks(1:PK.vector_number_complete_peaks(i),i)=PK.vector_left_normalized_amplitudes_complete_peaks;
            PK.matrix_right_amplitudes_complete_peaks(1:PK.vector_number_complete_peaks(i),i)=PK.vector_right_amplitudes_complete_peaks; % in AU
            PK.matrix_right_normalized_amplitudes_complete_peaks(1:PK.vector_number_complete_peaks(i),i)=PK.vector_right_normalized_amplitudes_complete_peaks;
            PK.matrix_left_velocities_complete_peaks(1:PK.vector_number_complete_peaks(i),i)=PK.vector_left_velocities_complete_peaks; % in AU/s
            PK.matrix_right_velocities_complete_peaks(1:PK.vector_number_complete_peaks(i),i)=PK.vector_right_velocities_complete_peaks; % in AU/s
            PK.vector_frequency_complete_peaks(i)=PK.vector_number_complete_peaks(i)/max(PK.vector_time); % in pk/s
            PK.vector_frequency_complete_valleys(i)=PK.vector_number_complete_valleys(i)/max(PK.vector_time); % in pk/s
            
            PK.vector_total_areas_complete_peaks=zeros(1,PK.vector_number_complete_peaks(i)); % in AU*s
            PK.vector_effective_areas_complete_peaks=zeros(1,PK.vector_number_complete_peaks(i)); % in AU*
            for m=1:PK.vector_number_complete_peaks(i)
                number_start_peak=find(PK.vector_time == PK.vector_time_complete_valleys(m));
                number_end_peak=find(PK.vector_time == PK.vector_time_complete_valleys(m+1));
                position=number_start_peak;
                local_area=0; % in AU*s
                while position<number_end_peak
                    local_area=(PK.matrix_filtered_fluorescences(position+1,i)+PK.matrix_filtered_fluorescences(position,i))*(PK.vector_time(position+1)-PK.vector_time(position))/2; % in AU*s
                    PK.vector_total_areas_complete_peaks(m)=PK.vector_total_areas_complete_peaks(m)+local_area; % in AU*s
                    position=position+1;
                end
                PK.vector_effective_areas_complete_peaks(m)=PK.vector_total_areas_complete_peaks(m)-(PK.matrix_filtered_fluorescences(number_end_peak,i)+PK.matrix_filtered_fluorescences(number_start_peak,i))*(PK.vector_time(number_end_peak)-PK.vector_time(number_start_peak))/2;
            end
            PK.matrix_effective_areas_complete_peaks(1:PK.vector_number_complete_peaks(i),i)=PK.vector_effective_areas_complete_peaks; % in AU*s
            
            PK.vector_maximum_left_amplitude_complete_peaks(i)=max(PK.vector_left_amplitudes_complete_peaks); % in AU
            PK.vector_maximum_left_normalized_amplitude_complete_peaks(i)=max(PK.vector_left_normalized_amplitudes_complete_peaks);
            PK.vector_maximum_right_amplitude_complete_peaks(i)=max(PK.vector_right_amplitudes_complete_peaks); % in AU
            PK.vector_maximum_right_normalized_amplitude_complete_peaks(i)=max(PK.vector_right_normalized_amplitudes_complete_peaks);
            %             k=1;
            %             while k==1
            %                 PK.vector_threshold_small_medium_peaks(i)=str2num(input('Enter the threshold of differentiation between small and medium peaks (between 0 and 1): ', 's')); % in AU
            %                 if PK.vector_threshold_small_medium_peaks(i)<=0 || PK.vector_threshold_small_medium_peaks(i)>=1
            %                     disp('The threshold of differentiation between small and medium peaks should be between 0 and 1');
            %                     n=2;
            %                 else
            %                     PK.vector_threshold_medium_large_peaks(i)=str2num(input('Enter the threshold of differentiation between medium and large peaks (between 0 and 1): ', 's')); % in AU
            %                     if PK.vector_threshold_medium_large_peaks(i)<=0 || PK.vector_threshold_medium_large_peaks(i)>=1 || PK.vector_threshold_small_medium_peaks(i)>=PK.vector_threshold_medium_large_peaks(i)
            %                         disp('The threshold of differentiation between medium and large peaks should be between 0 and 1 and greater than that between small and medium peaks');
            %                     else
            %                         k=2;
            %                     end
            %                 end
            %             end
            
            PK.vector_numbers_small_complete_peaks=find(PK.vector_left_amplitudes_complete_peaks<(PK.vector_threshold_small_medium_peaks(i)*PK.vector_maximum_left_amplitude_complete_peaks(i)));
            PK.vector_number_small_complete_peaks(i)=length(PK.vector_numbers_small_complete_peaks);
            PK.vector_left_durations_small_complete_peaks=PK.vector_left_durations_complete_peaks(PK.vector_numbers_small_complete_peaks); % in s
            PK.vector_right_durations_small_complete_peaks=PK.vector_right_durations_complete_peaks(PK.vector_numbers_small_complete_peaks); % in s
            PK.vector_left_amplitudes_small_complete_peaks=PK.vector_left_amplitudes_complete_peaks(PK.vector_numbers_small_complete_peaks); % in AU
            PK.vector_left_normalized_amplitudes_small_complete_peaks=PK.vector_left_normalized_amplitudes_complete_peaks(PK.vector_numbers_small_complete_peaks);
            PK.vector_right_amplitudes_small_complete_peaks=PK.vector_right_amplitudes_complete_peaks(PK.vector_numbers_small_complete_peaks); % in AU
            PK.vector_right_normalized_amplitudes_small_complete_peaks=PK.vector_right_normalized_amplitudes_complete_peaks(PK.vector_numbers_small_complete_peaks);
            PK.vector_left_velocities_small_complete_peaks=PK.vector_left_velocities_complete_peaks(PK.vector_numbers_small_complete_peaks); % in AU/s
            PK.vector_right_velocities_small_complete_peaks=PK.vector_right_velocities_complete_peaks(PK.vector_numbers_small_complete_peaks); % in AU/s
            PK.vector_areas_small_complete_peaks=PK.vector_effective_areas_complete_peaks(PK.vector_numbers_small_complete_peaks); % in AU*s
            PK.vector_frequency_small_complete_peaks(i)=PK.vector_number_small_complete_peaks(i)/max(PK.vector_time); % in pk/s
            PK.matrix_numbers_small_complete_peaks(1:PK.vector_number_small_complete_peaks(i),i)=PK.vector_numbers_small_complete_peaks;
            PK.matrix_left_durations_small_complete_peaks(1:PK.vector_number_small_complete_peaks(i),i)=PK.vector_left_durations_small_complete_peaks; % in s
            PK.matrix_right_durations_small_complete_peaks(1:PK.vector_number_small_complete_peaks(i),i)=PK.vector_right_durations_small_complete_peaks; % in s
            PK.matrix_left_amplitudes_small_complete_peaks(1:PK.vector_number_small_complete_peaks(i),i)=PK.vector_left_amplitudes_small_complete_peaks; % in AU
            PK.matrix_left_normalized_amplitudes_small_complete_peaks(1:PK.vector_number_small_complete_peaks(i),i)=PK.vector_left_normalized_amplitudes_small_complete_peaks;
            PK.matrix_right_amplitudes_small_complete_peaks(1:PK.vector_number_small_complete_peaks(i),i)=PK.vector_right_amplitudes_small_complete_peaks; % in AU
            PK.matrix_right_normalized_amplitudes_small_complete_peaks(1:PK.vector_number_small_complete_peaks(i),i)=PK.vector_right_normalized_amplitudes_small_complete_peaks;
            PK.matrix_left_velocities_small_complete_peaks(1:PK.vector_number_small_complete_peaks(i),i)=PK.vector_left_velocities_small_complete_peaks; % in AU/s
            PK.matrix_right_velocities_small_complete_peaks(1:PK.vector_number_small_complete_peaks(i),i)=PK.vector_right_velocities_small_complete_peaks; % in AU/s
            PK.matrix_areas_small_complete_peaks(1:PK.vector_number_small_complete_peaks(i),i)=PK.vector_areas_small_complete_peaks; % in AU*s
            
            PK.vector_numbers_medium_complete_peaks=find(PK.vector_left_amplitudes_complete_peaks>=(PK.vector_threshold_small_medium_peaks(i)*PK.vector_maximum_left_amplitude_complete_peaks(i)) & PK.vector_left_amplitudes_complete_peaks<(PK.vector_threshold_medium_large_peaks(i)*PK.vector_maximum_left_amplitude_complete_peaks(i)));
            PK.vector_number_medium_complete_peaks(i)=length(PK.vector_numbers_medium_complete_peaks);
            PK.vector_left_durations_medium_complete_peaks=PK.vector_left_durations_complete_peaks(PK.vector_numbers_medium_complete_peaks); % in s
            PK.vector_right_durations_medium_complete_peaks=PK.vector_right_durations_complete_peaks(PK.vector_numbers_medium_complete_peaks); % in s
            PK.vector_left_amplitudes_medium_complete_peaks=PK.vector_left_amplitudes_complete_peaks(PK.vector_numbers_medium_complete_peaks); % in AU
            PK.vector_left_normalized_amplitudes_medium_complete_peaks=PK.vector_left_normalized_amplitudes_complete_peaks(PK.vector_numbers_medium_complete_peaks);
            PK.vector_right_amplitudes_medium_complete_peaks=PK.vector_right_amplitudes_complete_peaks(PK.vector_numbers_medium_complete_peaks); % in AU
            PK.vector_right_normalized_amplitudes_medium_complete_peaks=PK.vector_right_normalized_amplitudes_complete_peaks(PK.vector_numbers_medium_complete_peaks);
            PK.vector_left_velocities_medium_complete_peaks=PK.vector_left_velocities_complete_peaks(PK.vector_numbers_medium_complete_peaks); % in AU/s
            PK.vector_right_velocities_medium_complete_peaks=PK.vector_right_velocities_complete_peaks(PK.vector_numbers_medium_complete_peaks); % in AU/s
            PK.vector_areas_medium_complete_peaks=PK.vector_effective_areas_complete_peaks(PK.vector_numbers_medium_complete_peaks); % in AU*s
            PK.vector_frequency_medium_complete_peaks(i)=PK.vector_number_medium_complete_peaks(i)/max(PK.vector_time); % in pk/s
            PK.matrix_numbers_medium_complete_peaks(1:PK.vector_number_medium_complete_peaks(i),i)=PK.vector_numbers_medium_complete_peaks;
            PK.matrix_left_durations_medium_complete_peaks(1:PK.vector_number_medium_complete_peaks(i),i)=PK.vector_left_durations_medium_complete_peaks; % in s
            PK.matrix_right_durations_medium_complete_peaks(1:PK.vector_number_medium_complete_peaks(i),i)=PK.vector_right_durations_medium_complete_peaks; % in s
            PK.matrix_left_amplitudes_medium_complete_peaks(1:PK.vector_number_medium_complete_peaks(i),i)=PK.vector_left_amplitudes_medium_complete_peaks; % in AU
            PK.matrix_left_normalized_amplitudes_medium_complete_peaks(1:PK.vector_number_medium_complete_peaks(i),i)=PK.vector_left_normalized_amplitudes_medium_complete_peaks;
            PK.matrix_right_amplitudes_medium_complete_peaks(1:PK.vector_number_medium_complete_peaks(i),i)=PK.vector_right_amplitudes_medium_complete_peaks; % in AU
            PK.matrix_right_normalized_amplitudes_medium_complete_peaks(1:PK.vector_number_medium_complete_peaks(i),i)=PK.vector_right_normalized_amplitudes_medium_complete_peaks;
            PK.matrix_left_velocities_medium_complete_peaks(1:PK.vector_number_medium_complete_peaks(i),i)=PK.vector_left_velocities_medium_complete_peaks; % in AU/s
            PK.matrix_right_velocities_medium_complete_peaks(1:PK.vector_number_medium_complete_peaks(i),i)=PK.vector_right_velocities_medium_complete_peaks; % in AU/s
            PK.matrix_areas_medium_complete_peaks(1:PK.vector_number_medium_complete_peaks(i),i)=PK.vector_areas_medium_complete_peaks; % in AU*s
            
            PK.vector_numbers_large_complete_peaks=find(PK.vector_left_amplitudes_complete_peaks>=(PK.vector_threshold_medium_large_peaks(i)*PK.vector_maximum_left_amplitude_complete_peaks(i)));
            PK.vector_number_large_complete_peaks(i)=length(PK.vector_numbers_large_complete_peaks);
            PK.vector_left_durations_large_complete_peaks=PK.vector_left_durations_complete_peaks(PK.vector_numbers_large_complete_peaks); % in s
            PK.vector_right_durations_large_complete_peaks=PK.vector_right_durations_complete_peaks(PK.vector_numbers_large_complete_peaks); % in s
            PK.vector_left_amplitudes_large_complete_peaks=PK.vector_left_amplitudes_complete_peaks(PK.vector_numbers_large_complete_peaks); % in AU
            PK.vector_left_normalized_amplitudes_large_complete_peaks=PK.vector_left_normalized_amplitudes_complete_peaks(PK.vector_numbers_large_complete_peaks);
            PK.vector_right_amplitudes_large_complete_peaks=PK.vector_right_amplitudes_complete_peaks(PK.vector_numbers_large_complete_peaks); % in AU
            PK.vector_right_normalized_amplitudes_large_complete_peaks=PK.vector_right_normalized_amplitudes_complete_peaks(PK.vector_numbers_large_complete_peaks);
            PK.vector_left_velocities_large_complete_peaks=PK.vector_left_velocities_complete_peaks(PK.vector_numbers_large_complete_peaks); % in AU/s
            PK.vector_right_velocities_large_complete_peaks=PK.vector_right_velocities_complete_peaks(PK.vector_numbers_large_complete_peaks); % in AU/s
            PK.vector_areas_large_complete_peaks=PK.vector_effective_areas_complete_peaks(PK.vector_numbers_large_complete_peaks); % in AU*s
            PK.vector_frequency_large_complete_peaks(i)=PK.vector_number_large_complete_peaks(i)/max(PK.vector_time); % in pk/s
            PK.matrix_numbers_large_complete_peaks(1:PK.vector_number_large_complete_peaks(i),i)=PK.vector_numbers_large_complete_peaks;
            PK.matrix_left_durations_large_complete_peaks(1:PK.vector_number_large_complete_peaks(i),i)=PK.vector_left_durations_large_complete_peaks; % in s
            PK.matrix_right_durations_large_complete_peaks(1:PK.vector_number_large_complete_peaks(i),i)=PK.vector_right_durations_large_complete_peaks; % in s
            PK.matrix_left_amplitudes_large_complete_peaks(1:PK.vector_number_large_complete_peaks(i),i)=PK.vector_left_amplitudes_large_complete_peaks; % in AU
            PK.matrix_left_normalized_amplitudes_large_complete_peaks(1:PK.vector_number_large_complete_peaks(i),i)=PK.vector_left_normalized_amplitudes_large_complete_peaks;
            PK.matrix_right_amplitudes_large_complete_peaks(1:PK.vector_number_large_complete_peaks(i),i)=PK.vector_right_amplitudes_large_complete_peaks; % in AU
            PK.matrix_right_normalized_amplitudes_large_complete_peaks(1:PK.vector_number_large_complete_peaks(i),i)=PK.vector_right_normalized_amplitudes_large_complete_peaks;
            PK.matrix_left_velocities_large_complete_peaks(1:PK.vector_number_large_complete_peaks(i),i)=PK.vector_left_velocities_large_complete_peaks; % in AU/s
            PK.matrix_right_velocities_large_complete_peaks(1:PK.vector_number_large_complete_peaks(i),i)=PK.vector_right_velocities_large_complete_peaks; % in AU/s
            PK.matrix_areas_large_complete_peaks(1:PK.vector_number_large_complete_peaks(i),i)=PK.vector_areas_large_complete_peaks; % in AU*s
            
            
        end
        
        function PK = CalculateParameters(PK,varargin)
            p=inputParser;
            addRequired(p,'i')
            parse(p, varargin{1:end});
            i=p.Results.i;
            dper=[0.8 0.7 0.5 0.3 0.1];
            PK.posm(:,i)=nan*ones(100,1);
            PK.posM(:,i)=nan*ones(100,1);
            PK.ms(:,i)=nan*ones(100,1);
            PK.Ms(:,i)=nan*ones(100,1);
            PK.xMr(:,i)=nan*ones(100,1);
            PK.xmc(:,i)=nan*ones(100,1);
            PK.xMc(:,i)=nan*ones(100,1);
            PK.mmvg(:,i)=nan*ones(100,1);
            PK.mmvd(:,i)=nan*ones(100,1);
            PK.M(:,i)=nan*ones(100,1);
            PK.Tauc(:,i)=nan*ones(100,1);
            PK.Taur(:,i)=nan*ones(100,1);
            PK.Taud(:,i)=nan*ones(100,1);
            PK.Taus(:,i)=nan*ones(100,1);
            PK.Tauc(:,i)=nan*ones(100,1);
            PK.pc(:,i)=nan*ones(100,1);
            PK.pr(:,i)=nan*ones(100,1);
            PK.hc(:,i)=nan*ones(100,1);
            PK.hr(:,i)=nan*ones(100,1);
            PK.Aire(:,i)=nan*ones(100,1);
            PK.locc(:,i)=nan*ones(100,1);
%             PK.locr(:,i)=nan*ones(100,1);

            PK.Mper(:,:,i)=nan*ones(100,5,1);
            PK.posper(:,:,i)=nan*ones(100,5,1);
            
            
            
            
            
            win=10;
            sl=PK.sm(i);
            Signal=PK.matrix_filtered_fluorescences(~isnan(PK.matrix_filtered_fluorescences(:,i)),i);
            %             Signal=    PK.matrix_rough_fluorescences(:,i);
            
            
            ss=smooth(Signal,sl,'loess');
            
            
            
            dd2=gradient(ss);
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
                            %                                                                                 plot(locc(ii),PK.Area(kk).FinalSignal(locc(ii)),'+r')
                            %                                                                                 plot(locr(pl-1),PK.Area(kk).FinalSignal(locr(pl-1)),'+g')
                            %                                 %
                            [M(count),ll(ii)]=max(Signal(locc(ii):locr(pl)));
                            ll(ii)=ll(ii)+locc(ii)-1;
                            
                            % [mmvg(count),vvg(ii)]=min(PK.Area(kk).FinalSignal(locr(pl-1):locc(ii)));
                            [mmvg(count),vvg(ii)]=min(ss(locr(pl-1):locc(ii)));
                            vvg(ii)=vvg(ii)+locr(pl-1)-1;
                            %             [mmvd(count),vvd(ii)]=min(PK.Area(kk).FinalSignal(locr(pl):locc(ii+1)));
                            [mmvd(count),vvd(ii)]=min(ss(locr(pl):locc(ii+1)));
                            vvd(ii)=vvd(ii)+locr(pl)-1;
                            %                                                                                plot(ll(ii),M(ii),'+b')
                            %                                                                                 plot(vv(ii),mmv(ii),'+y')
                            
                            
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
                            PK.pc(count,i)=maxc(ii)*PK.PixelSize*PK.SamplingFrequency;
                            PK.pr(count,i)=minr(pl)*PK.PixelSize*PK.SamplingFrequency;
                            PK.hr(count,i)=-minr(pl)*(xmr(count)-xMr(count))*PK.PixelSize;
                            PK.hc(count,i)=maxc(ii)*(xMc(count)-xmc(count))*PK.PixelSize;
                            PK.M(count,i)=M(count)*PK.PixelSize;
                            PK.mmvg(count,i)=mmvg(count)*PK.PixelSize;
                            PK.mmvd(count,i)=mmvd(count)*PK.PixelSize;
                            [Ms(count),posM(count)]=max(Signal(xmc(count):xmr(count)));
                            posM(count)=posM(count)+xmc(count)-1;
                            
                            count=count+1;
                        end
                    end
                end
            end
            for kk=1:length(xMr)-1
                [ms(kk),posm(kk)]=min(Signal(xMr(kk):xMr(kk+1)));
                posm(kk)=posm(kk)+xMr(kk)-1;
               
            end
            [ms(length(xMr)),posm(length(xMr))]=min(Signal(xMr(length(xMr)):end));
            posm(length(xMr))=posm(length(xMr))+xMr(length(xMr))-1;
         
            for uu=1:5
            for kk=1:length(xMr(:))
                dd=dper(uu)*(Ms(kk)-ms(kk));
                posM(kk)
                tutu=find(diff(Signal(posM(kk):posm(kk))-ms(kk)>dd)==-1);
                if ~isempty(tutu)
                posper(kk,uu)=tutu(1);
                posper(kk,uu)=posper(kk,uu)+posM(kk)-1;
                Mper(kk,uu)=ms(kk)+dd;
               
                end
            end
            end
            
            PK.posper(1:size(posper,1),:,i)=posper/PK.SamplingFrequency;
            PK.Mper(1:size(posper,1),:,i)=Mper*PK.PixelSize;
            PK.posm(1:length(xmc),i)=posm/PK.SamplingFrequency;
            PK.posM(1:length(xmc),i)=posM/PK.SamplingFrequency;
            PK.ms(1:length(xmc),i)=ms*PK.PixelSize;
            PK.Ms(1:length(xmc),i)=Ms*PK.PixelSize;
            PK.xmc(1:length(xmc),i)=xmc/PK.SamplingFrequency;
            PK.xMc(1:length(xMc),i)=xMc/PK.SamplingFrequency;
            PK.xmr(1:length(xmr),i)=xmr/PK.SamplingFrequency;
            PK.xMr(1:length(xMr),i)=xMr/PK.SamplingFrequency;
            PK.locc(1:length(locc),i)=locc/PK.SamplingFrequency;
            PK.locr(1:length(locr),i)=locr/PK.SamplingFrequency;
            
            
%             PK.MedTauc = median(PK.Tauc);
%             PK.StdTauc = std(PK.Tauc);
%             PK.MedTaur = median(PK.Taur);
%             PK.StdTaur = std(PK.Taur);
%             PK.MedTaus = median(PK.Taus);
%             PK.StdTaus = std(PK.Taus);
%             PK.MedTausBaz = median(PK.TausBaz);
%             PK.StdTausBaz = std(PK.TausBaz);
%             PK.MedTaud = median(PK.Taud);
%             PK.StdTaud = std(PK.Taud);
%             PK.MedTaudBaz = median(PK.TaudBaz);
%             PK.StdTaudBaz = std(PK.TaudBaz);
%             PK.MedAire = median(PK.Aire);
%             PK.StdAire = std(PK.Aire);
%             PK.Medpc = median(PK.pc);
%             PK.Stdpc = std(PK.pc);
%             PK.Medpr = median(PK.pr);
%             PK.Stdpr = std(PK.pr);
%             PK.Medhr = median(PK.hr);
%             PK.Stdhr = std(PK.hr);
%             PK.Medhc = median(PK.hc);
%             PK.Stdhc = std(PK.hc);
%             PK.MedM = median(PK.M);
%             PK.StdM = std(PK.M);
%             PK.Medmmvg = median(PK.mmvg);
%             PK.Stdmmvg = std(PK.mmvg);
%             PK.Medmmvd = median(PK.mmvd);
%             PK.Stdmmvd = std(PK.mmvd);
        end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        function PK=Save(PK,varargin)
            results_pathname=varargin{1};
            
            
            
            for x=1:PK.number_cells
                
                clear Tf
                Tf= array2table(NaN*ones(100,16));
                Tf.Properties.VariableNames = {'peak_time', 'freq_small_Ca_events',...
                    'freq_medium_Ca_events','frequency_of_Ca_transients','ascending_time','decay_time','decay_time_90','decay_time_70','decay_time_50','decay_time_30','decay_time_20',...
                    'abs_amp_of_Ca_trans','abs_amp_max','norm_amp_of_Ca_trans','norm_amp_max','velocity'};
                
                
                
                val=PK.matrix_time_complete_peaks(PK.matrix_time_complete_peaks(:,x)>0,x);
                Tf.peak_time(1:length(val))=val;
                val=PK.vector_frequency_small_complete_peaks(:,x);
                Tf.freq_small_Ca_events(1:length(val))=val;
                val=PK.vector_frequency_medium_complete_peaks(:,x);
                Tf.freq_medium_Ca_events(1:length(val))=val;
                val=PK.vector_frequency_large_complete_peaks(:,x);
                Tf.frequency_of_Ca_transients(1:length(val))=val;
                val=PK.matrix_left_durations_large_complete_peaks(PK.matrix_left_durations_large_complete_peaks(:,x)>0,x);
                Tf.ascending_time(1:length(val))=val;
                val=PK.matrix_right_durations_large_complete_peaks(PK.matrix_right_durations_large_complete_peaks(:,x)>0,x);
                Tf.decay_time(1:length(val))=val;
                val=0.9*Tf.decay_time;
                Tf.decay_time_90(1:length(val))=val;
                val=0.7*Tf.decay_time;
                Tf.decay_time_70(1:length(val))=val;
                val=0.5*Tf.decay_time;
                Tf.decay_time_50(1:length(val))=val;
                val=0.2*Tf.decay_time;
                Tf.decay_time_20(1:length(val))=val;
                val=PK.matrix_left_amplitudes_large_complete_peaks(PK.matrix_left_amplitudes_large_complete_peaks(:,x)>0,x);
                Tf.abs_amp_of_Ca_trans(1:length(val))=val;
                val=0.3*Tf.decay_time;
                Tf.decay_time_30(1:length(val))=val;
                val=PK.vector_maximum_left_amplitude_complete_peaks(:,x);
                Tf.abs_amp_max(1:length(val))=val;
                val=PK.matrix_left_normalized_amplitudes_large_complete_peaks(PK.matrix_left_normalized_amplitudes_large_complete_peaks(:,x)>0,x);
                Tf.norm_amp_of_Ca_trans(1:length(val))=val;
                val=PK.vector_maximum_left_normalized_amplitude_complete_peaks(:,x);
                Tf.norm_amp_max(1:length(val))=val;
                val=PK.matrix_left_velocities_large_complete_peaks(PK.matrix_left_velocities_large_complete_peaks(:,x)>0,x);
                Tf.velocity(1:length(val))=val;
                
                writetable(Tf,results_pathname,'sheet',['cell' num2str(x)])
                
            end
            
        end
        
        function PK=SaveOne(PK,varargin)
            
        end
        
    end
    
end