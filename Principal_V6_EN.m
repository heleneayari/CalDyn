clear
close all;

operating_system=menu('Which operating system do you use?','Macintosh','Windows');
if operating_system==1
    rough_data_foldername='Rough_data/';
else
    rough_data_foldername='Rough_data\';
end 
rough_data_filename=input('Enter the rough data filename (without extension): ', 's');
rough_data_pathname=[rough_data_foldername rough_data_filename];
matrix_rough_data=xlsread([rough_data_pathname '.xlsx']);
number_acquisitions=size(matrix_rough_data,1);
number_cells=size(matrix_rough_data,2)-1;
vector_time=matrix_rough_data(:,1);
matrix_rough_fluorescences=matrix_rough_data(:,2:(number_cells+1));
vector_filtering_polynomial_order=zeros(1,number_cells); 
vector_filtering_frame_length=zeros(1,number_cells); 
matrix_filtered_fluorescences=zeros(number_acquisitions,number_cells);
vector_threshold_peaks_detection=zeros(1,number_cells);
vector_number_peaks=zeros(1,number_cells);
vector_number_valleys=zeros(1,number_cells);
matrix_time_peaks=zeros(number_acquisitions,number_cells);
matrix_time_valleys=zeros(number_acquisitions,number_cells);
matrix_fluorescences_peaks=zeros(number_acquisitions,number_cells);
matrix_fluorescences_valleys=zeros(number_acquisitions,number_cells);
vector_number_complete_peaks=zeros(1,number_cells);
vector_number_complete_valleys=zeros(1,number_cells);
matrix_time_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_time_complete_valleys=zeros(number_acquisitions,number_cells);
matrix_fluorescences_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_fluorescences_complete_valleys=zeros(number_acquisitions,number_cells);
matrix_left_durations_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_right_durations_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_left_amplitudes_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_left_normalized_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_right_amplitudes_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_right_normalized_amplitudes_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_left_velocities_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_right_velocities_complete_peaks=zeros(number_acquisitions,number_cells);
vector_frequency_complete_peaks=zeros(1,number_cells);
vector_frequency_complete_valleys=zeros(1,number_cells);
matrix_effective_areas_complete_peaks=zeros(number_acquisitions,number_cells);
vector_maximum_left_amplitude_complete_peaks=zeros(1,number_cells);
vector_maximum_left_normalized_amplitude_complete_peaks=zeros(1,number_cells);
vector_maximum_right_amplitude_complete_peaks=zeros(1,number_cells);
vector_maximum_right_normalized_amplitude_complete_peaks=zeros(1,number_cells);
vector_threshold_small_medium_peaks=zeros(1,number_cells);
vector_threshold_medium_large_peaks=zeros(1,number_cells);
vector_number_small_complete_peaks=zeros(1,number_cells);
vector_frequency_small_complete_peaks=zeros(1,number_cells);
matrix_numbers_small_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_left_durations_small_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_right_durations_small_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_left_amplitudes_small_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_left_normalized_amplitudes_small_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_right_amplitudes_small_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_right_normalized_amplitudes_small_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_left_velocities_small_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_right_velocities_small_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_areas_small_complete_peaks=zeros(number_acquisitions,number_cells);
vector_number_medium_complete_peaks=zeros(1,number_cells);
vector_frequency_medium_complete_peaks=zeros(1,number_cells);
matrix_numbers_medium_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_left_durations_medium_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_right_durations_medium_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_left_amplitudes_medium_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_left_normalized_amplitudes_medium_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_right_amplitudes_medium_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_right_normalized_amplitudes_medium_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_left_velocities_medium_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_right_velocities_medium_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_areas_medium_complete_peaks=zeros(number_acquisitions,number_cells);
vector_number_large_complete_peaks=zeros(1,number_cells); %%%%%%
vector_frequency_large_complete_peaks=zeros(1,number_cells); %%%%%%
matrix_numbers_large_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_left_durations_large_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_right_durations_large_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_left_amplitudes_large_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_left_normalized_amplitudes_large_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_right_amplitudes_large_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_right_normalized_amplitudes_large_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_left_velocities_large_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_right_velocities_large_complete_peaks=zeros(number_acquisitions,number_cells);
matrix_areas_large_complete_peaks=zeros(number_acquisitions,number_cells);

for i=1:number_cells
    j=1;
    s=1;
    while j==1
        vector_filtering_polynomial_order(i)=str2num(input(['Enter the polynomial order of Savitzky-Golay filtering for cell ' num2str(i) ': '], 's')); 
        vector_filtering_frame_length(i)=str2num(input(['Enter the frame length of Savitzky-Golay filtering for cell ' num2str(i) ': '], 's')); 
        matrix_filtered_fluorescences(:,i)=sgolayfilt(matrix_rough_fluorescences(:,i),vector_filtering_polynomial_order(i),vector_filtering_frame_length(i)); % in AU
        
        close(figure(1));
        figure(1)
        subplot(3,1,1);
        plot(vector_time,matrix_rough_fluorescences(:,i),'-b');
        grid on
        title('(a) Rough signal','interpreter','latex');
        xlabel('Time (s)','interpreter','latex');
        ylabel('Fluorescence (AU)','interpreter','latex');
        legend('Experimental data (rough)');
        subplot(3,1,2);
        plot(vector_time,matrix_filtered_fluorescences(:,i),'-k');
        grid on
        title('(b) Filtered signal','interpreter','latex')
        xlabel('Time (s)','interpreter','latex');
        ylabel('Fluorescence (AU)','interpreter','latex');
        legend('Experimental data (filtered)');
        subplot(3,1,3);
        plot(vector_time,matrix_rough_fluorescences(:,i),'-b');
        hold on
        plot(vector_time,matrix_filtered_fluorescences(:,i),'-k');
        grid on
        title('(c) Superposition of rough and filtered signals','interpreter','latex');
        xlabel('Time (s)','interpreter','latex');
        ylabel('Fluorescence (AU)','interpreter','latex');
        legend('Experimental data (rough)','Experimental data (filtered)');
        
        j=menu('Do you want to','change the polynomial order and/or the frame length?','continue the data processing?');
        if j==2
            while s==1
                vector_threshold_peaks_detection(i)=input('Look to Figure 1b and enter the peak detection threshold (minimum amplitude): '); % in AU
                [vector_fluorescences_peaks,vector_time_peaks]=findpeaks(matrix_filtered_fluorescences(:,i),vector_time,'MinPeakProminence',vector_threshold_peaks_detection(i));
                [vector_fluorescences_valleys,vector_time_valleys]=findpeaks(-matrix_filtered_fluorescences(:,i),vector_time,'MinPeakProminence',vector_threshold_peaks_detection(i));
                if length(vector_time_peaks)==1 && isempty(vector_time_valleys)
                    vector_fluorescences_valleys(1,1)=max(-matrix_filtered_fluorescences(1:find(vector_time==vector_time_peaks),1));
                    vector_time_valleys(1,1)=vector_time(-matrix_filtered_fluorescences(:,i)==vector_fluorescences_valleys(1,1));
                    vector_fluorescences_valleys(2,1)=max(-matrix_filtered_fluorescences(find(vector_time==vector_time_peaks):length(vector_time),1));
                    vector_time_valleys(2,1)=vector_time(-matrix_filtered_fluorescences(:,i)==vector_fluorescences_valleys(2,1));
                end
                vector_fluorescences_valleys=-vector_fluorescences_valleys; % in AU
                vector_number_peaks(i)=length(vector_time_peaks);
                vector_number_valleys(i)=length(vector_time_valleys);
                matrix_time_peaks(1:vector_number_peaks(i),i)=vector_time_peaks; % in s
                matrix_time_valleys(1:vector_number_valleys(i),i)=vector_time_valleys; % in s
                matrix_fluorescences_peaks(1:vector_number_peaks(i),i)=vector_fluorescences_peaks; % in AU 
                matrix_fluorescences_valleys(1:vector_number_valleys(i),i)=vector_fluorescences_valleys; % in AU
                
                close(figure(2));
                figure(2)
                plot(vector_time,matrix_filtered_fluorescences(:,i),'-k');
                grid on
                title('Filtered signal','interpreter','latex')
                xlabel('Time (s)','interpreter','latex');
                ylabel('Fluorescence (AU)','interpreter','latex');
                hold on
                plot(vector_time_peaks,vector_fluorescences_peaks,'+r');
                hold on
                plot(vector_time_valleys,vector_fluorescences_valleys,'xb');
                legend('Experimental data (filtered)','Peaks','Valleys');
                s=menu('Do you want to','change the peak detection threshold?','continue the data processing?');
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Caracterisation des peaks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if vector_time_valleys(1)<vector_time_peaks(1)
                if vector_number_peaks(i)<vector_number_valleys(i)
                    vector_number_complete_peaks(i)=vector_number_peaks(i);
                    vector_number_complete_valleys(i)=vector_number_valleys(i);
                    vector_time_complete_peaks=vector_time_peaks; % in s
                    vector_time_complete_valleys=vector_time_valleys; % in s
                    vector_fluorescences_complete_peaks=vector_fluorescences_peaks; % in AU
                    vector_fluorescences_complete_valleys=vector_fluorescences_valleys; % in AU
                else
                    vector_number_complete_peaks(i)=vector_number_peaks(i)-1;
                    vector_number_complete_valleys(i)=vector_number_valleys(i);
                    vector_time_complete_peaks=vector_time_peaks(1:vector_number_complete_peaks(i)); % in s
                    vector_time_complete_valleys=vector_time_valleys; % in s
                    vector_fluorescences_complete_peaks=vector_fluorescences_peaks(1:vector_number_complete_peaks(i)); % in AU
                    vector_fluorescences_complete_valleys=vector_fluorescences_valleys; % in AU
                end
            else
                if vector_number_peaks(i)==vector_number_valleys(i)
                    vector_number_complete_peaks(i)=vector_number_peaks(i)-1;
                    vector_number_complete_valleys(i)=vector_number_valleys(i);
                    vector_time_complete_peaks=vector_time_peaks(2:(vector_number_complete_peaks(i)+1)); % in s
                    vector_time_complete_valleys=vector_time_valleys; % in s
                    vector_fluorescences_complete_peaks=vector_fluorescences_peaks(2:(vector_number_complete_peaks(i)+1)); % in AU
                    vector_fluorescences_complete_valleys=vector_fluorescences_valleys; % in AU
                else
                    vector_number_complete_peaks(i)=vector_number_peaks(i)-2;
                    vector_number_complete_valleys(i)=vector_number_valleys(i);
                    vector_time_complete_peaks=vector_time_peaks(2:(vector_number_complete_peaks(i)+1)); % in s
                    vector_time_complete_valleys=vector_time_valleys; % in s
                    vector_fluorescences_complete_peaks=vector_fluorescences_peaks(2:(vector_number_complete_peaks(i)+1)); % in AU
                    vector_fluorescences_complete_valleys=vector_fluorescences_valleys; % in AU
                end
            end
            vector_left_durations_complete_peaks=vector_time_complete_peaks-vector_time_complete_valleys(1:(vector_number_complete_valleys(i)-1)); % in s
            vector_right_durations_complete_peaks=-(vector_time_complete_peaks-vector_time_complete_valleys(2:vector_number_complete_valleys(i))); % in s    
            vector_left_amplitudes_complete_peaks=vector_fluorescences_complete_peaks-vector_fluorescences_complete_valleys(1:(vector_number_complete_valleys(i)-1)); % in AU
            vector_left_normalized_amplitudes_complete_peaks=vector_left_amplitudes_complete_peaks./vector_fluorescences_complete_valleys(1:(vector_number_complete_valleys(i)-1));
            vector_right_amplitudes_complete_peaks=vector_fluorescences_complete_peaks-vector_fluorescences_complete_valleys(2:vector_number_complete_valleys(i)); % in AU
            vector_right_normalized_amplitudes_complete_peaks=vector_right_amplitudes_complete_peaks./vector_fluorescences_complete_valleys(2:vector_number_complete_valleys(i));
            vector_left_velocities_complete_peaks=vector_left_amplitudes_complete_peaks./vector_left_durations_complete_peaks; % in AU/s
            vector_right_velocities_complete_peaks=vector_right_amplitudes_complete_peaks./vector_right_durations_complete_peaks; % in AU/s
            matrix_time_complete_peaks(1:vector_number_complete_peaks(i),i)=vector_time_complete_peaks; % in s
            matrix_time_complete_valleys(1:vector_number_complete_valleys(i),i)=vector_time_complete_valleys; % in s
            matrix_fluorescences_complete_peaks(1:vector_number_complete_peaks(i),i)=vector_fluorescences_complete_peaks; % in AU
            matrix_fluorescences_complete_valleys(1:vector_number_complete_valleys(i),i)=vector_fluorescences_complete_valleys; % in AU
            matrix_left_durations_complete_peaks(1:vector_number_complete_peaks(i),i)=vector_left_durations_complete_peaks; % in s
            matrix_right_durations_complete_peaks(1:vector_number_complete_peaks(i),i)=vector_right_durations_complete_peaks; % in s
            matrix_left_amplitudes_complete_peaks(1:vector_number_complete_peaks(i),i)=vector_left_amplitudes_complete_peaks; % in AU
            matrix_left_normalized_complete_peaks(1:vector_number_complete_peaks(i),i)=vector_left_normalized_amplitudes_complete_peaks;
            matrix_right_amplitudes_complete_peaks(1:vector_number_complete_peaks(i),i)=vector_right_amplitudes_complete_peaks; % in AU
            matrix_right_normalized_amplitudes_complete_peaks(1:vector_number_complete_peaks(i),i)=vector_right_normalized_amplitudes_complete_peaks;
            matrix_left_velocities_complete_peaks(1:vector_number_complete_peaks(i),i)=vector_left_velocities_complete_peaks; % in AU/s
            matrix_right_velocities_complete_peaks(1:vector_number_complete_peaks(i),i)=vector_right_velocities_complete_peaks; % in AU/s
            vector_frequency_complete_peaks(i)=vector_number_complete_peaks(i)/max(vector_time); % in pk/s
            vector_frequency_complete_valleys(i)=vector_number_complete_valleys(i)/max(vector_time); % in pk/s
            
            vector_total_areas_complete_peaks=zeros(1,vector_number_complete_peaks(i)); % in AU*s
            vector_effective_areas_complete_peaks=zeros(1,vector_number_complete_peaks(i)); % in AU*
            for m=1:vector_number_complete_peaks(i)
                number_start_peak=find(vector_time == vector_time_complete_valleys(m));
                number_end_peak=find(vector_time == vector_time_complete_valleys(m+1));
                position=number_start_peak;
                local_area=0; % in AU*s
                while position<number_end_peak
                    local_area=(matrix_filtered_fluorescences(position+1,i)+matrix_filtered_fluorescences(position,i))*(vector_time(position+1)-vector_time(position))/2; % in AU*s
                    vector_total_areas_complete_peaks(m)=vector_total_areas_complete_peaks(m)+local_area; % in AU*s
                    position=position+1;
                end
                vector_effective_areas_complete_peaks(m)=vector_total_areas_complete_peaks(m)-(matrix_filtered_fluorescences(number_end_peak,i)+matrix_filtered_fluorescences(number_start_peak,i))*(vector_time(number_end_peak)-vector_time(number_start_peak))/2;
            end
            matrix_effective_areas_complete_peaks(1:vector_number_complete_peaks(i),i)=vector_effective_areas_complete_peaks; % in AU*s

            vector_maximum_left_amplitude_complete_peaks(i)=max(vector_left_amplitudes_complete_peaks); % in AU
            vector_maximum_left_normalized_amplitude_complete_peaks(i)=max(vector_left_normalized_amplitudes_complete_peaks);
            vector_maximum_right_amplitude_complete_peaks(i)=max(vector_right_amplitudes_complete_peaks); % in AU
            vector_maximum_right_normalized_amplitude_complete_peaks(i)=max(vector_right_normalized_amplitudes_complete_peaks);
            k=1;
            while k==1
                vector_threshold_small_medium_peaks(i)=str2num(input('Enter the threshold of differentiation between small and medium peaks (between 0 and 1): ', 's')); % in AU
                if vector_threshold_small_medium_peaks(i)<=0 || vector_threshold_small_medium_peaks(i)>=1
                    disp('The threshold of differentiation between small and medium peaks should be between 0 and 1');
                    n=2;
                else
                    vector_threshold_medium_large_peaks(i)=str2num(input('Enter the threshold of differentiation between medium and large peaks (between 0 and 1): ', 's')); % in AU
                    if vector_threshold_medium_large_peaks(i)<=0 || vector_threshold_medium_large_peaks(i)>=1 || vector_threshold_small_medium_peaks(i)>=vector_threshold_medium_large_peaks(i)
                        disp('The threshold of differentiation between medium and large peaks should be between 0 and 1 and greater than that between small and medium peaks');
                    else
                        k=2;
                    end
                end
            end
            
            vector_numbers_small_complete_peaks=find(vector_left_amplitudes_complete_peaks<(vector_threshold_small_medium_peaks(i)*vector_maximum_left_amplitude_complete_peaks(i)));
            vector_number_small_complete_peaks(i)=length(vector_numbers_small_complete_peaks);
            vector_left_durations_small_complete_peaks=vector_left_durations_complete_peaks(vector_numbers_small_complete_peaks); % in s
            vector_right_durations_small_complete_peaks=vector_right_durations_complete_peaks(vector_numbers_small_complete_peaks); % in s
            vector_left_amplitudes_small_complete_peaks=vector_left_amplitudes_complete_peaks(vector_numbers_small_complete_peaks); % in AU
            vector_left_normalized_amplitudes_small_complete_peaks=vector_left_normalized_amplitudes_complete_peaks(vector_numbers_small_complete_peaks);
            vector_right_amplitudes_small_complete_peaks=vector_right_amplitudes_complete_peaks(vector_numbers_small_complete_peaks); % in AU
            vector_right_normalized_amplitudes_small_complete_peaks=vector_right_normalized_amplitudes_complete_peaks(vector_numbers_small_complete_peaks);
            vector_left_velocities_small_complete_peaks=vector_left_velocities_complete_peaks(vector_numbers_small_complete_peaks); % in AU/s
            vector_right_velocities_small_complete_peaks=vector_right_velocities_complete_peaks(vector_numbers_small_complete_peaks); % in AU/s
            vector_areas_small_complete_peaks=vector_effective_areas_complete_peaks(vector_numbers_small_complete_peaks); % in AU*s
            vector_frequency_small_complete_peaks(i)=vector_number_small_complete_peaks(i)/max(vector_time); % in pk/s
            matrix_numbers_small_complete_peaks(1:vector_number_small_complete_peaks(i),i)=vector_numbers_small_complete_peaks;
            matrix_left_durations_small_complete_peaks(1:vector_number_small_complete_peaks(i),i)=vector_left_durations_small_complete_peaks; % in s
            matrix_right_durations_small_complete_peaks(1:vector_number_small_complete_peaks(i),i)=vector_right_durations_small_complete_peaks; % in s
            matrix_left_amplitudes_small_complete_peaks(1:vector_number_small_complete_peaks(i),i)=vector_left_amplitudes_small_complete_peaks; % in AU
            matrix_left_normalized_amplitudes_small_complete_peaks(1:vector_number_small_complete_peaks(i),i)=vector_left_normalized_amplitudes_small_complete_peaks;
            matrix_right_amplitudes_small_complete_peaks(1:vector_number_small_complete_peaks(i),i)=vector_right_amplitudes_small_complete_peaks; % in AU
            matrix_right_normalized_amplitudes_small_complete_peaks(1:vector_number_small_complete_peaks(i),i)=vector_right_normalized_amplitudes_small_complete_peaks;
            matrix_left_velocities_small_complete_peaks(1:vector_number_small_complete_peaks(i),i)=vector_left_velocities_small_complete_peaks; % in AU/s
            matrix_right_velocities_small_complete_peaks(1:vector_number_small_complete_peaks(i),i)=vector_right_velocities_small_complete_peaks; % in AU/s
            matrix_areas_small_complete_peaks(1:vector_number_small_complete_peaks(i),i)=vector_areas_small_complete_peaks; % in AU*s
            
            vector_numbers_medium_complete_peaks=find(vector_left_amplitudes_complete_peaks>=(vector_threshold_small_medium_peaks(i)*vector_maximum_left_amplitude_complete_peaks(i)) & vector_left_amplitudes_complete_peaks<(vector_threshold_medium_large_peaks(i)*vector_maximum_left_amplitude_complete_peaks(i)));
            vector_number_medium_complete_peaks(i)=length(vector_numbers_medium_complete_peaks);
            vector_left_durations_medium_complete_peaks=vector_left_durations_complete_peaks(vector_numbers_medium_complete_peaks); % in s
            vector_right_durations_medium_complete_peaks=vector_right_durations_complete_peaks(vector_numbers_medium_complete_peaks); % in s
            vector_left_amplitudes_medium_complete_peaks=vector_left_amplitudes_complete_peaks(vector_numbers_medium_complete_peaks); % in AU
            vector_left_normalized_amplitudes_medium_complete_peaks=vector_left_normalized_amplitudes_complete_peaks(vector_numbers_medium_complete_peaks);
            vector_right_amplitudes_medium_complete_peaks=vector_right_amplitudes_complete_peaks(vector_numbers_medium_complete_peaks); % in AU
            vector_right_normalized_amplitudes_medium_complete_peaks=vector_right_normalized_amplitudes_complete_peaks(vector_numbers_medium_complete_peaks);
            vector_left_velocities_medium_complete_peaks=vector_left_velocities_complete_peaks(vector_numbers_medium_complete_peaks); % in AU/s
            vector_right_velocities_medium_complete_peaks=vector_right_velocities_complete_peaks(vector_numbers_medium_complete_peaks); % in AU/s
            vector_areas_medium_complete_peaks=vector_effective_areas_complete_peaks(vector_numbers_medium_complete_peaks); % in AU*s            
            vector_frequency_medium_complete_peaks(i)=vector_number_medium_complete_peaks(i)/max(vector_time); % in pk/s
            matrix_numbers_medium_complete_peaks(1:vector_number_medium_complete_peaks(i),i)=vector_numbers_medium_complete_peaks;
            matrix_left_durations_medium_complete_peaks(1:vector_number_medium_complete_peaks(i),i)=vector_left_durations_medium_complete_peaks; % in s
            matrix_right_durations_medium_complete_peaks(1:vector_number_medium_complete_peaks(i),i)=vector_right_durations_medium_complete_peaks; % in s
            matrix_left_amplitudes_medium_complete_peaks(1:vector_number_medium_complete_peaks(i),i)=vector_left_amplitudes_medium_complete_peaks; % in AU
            matrix_left_normalized_amplitudes_medium_complete_peaks(1:vector_number_medium_complete_peaks(i),i)=vector_left_normalized_amplitudes_medium_complete_peaks;            
            matrix_right_amplitudes_medium_complete_peaks(1:vector_number_medium_complete_peaks(i),i)=vector_right_amplitudes_medium_complete_peaks; % in AU
            matrix_right_normalized_amplitudes_medium_complete_peaks(1:vector_number_medium_complete_peaks(i),i)=vector_right_normalized_amplitudes_medium_complete_peaks;
            matrix_left_velocities_medium_complete_peaks(1:vector_number_medium_complete_peaks(i),i)=vector_left_velocities_medium_complete_peaks; % in AU/s
            matrix_right_velocities_medium_complete_peaks(1:vector_number_medium_complete_peaks(i),i)=vector_right_velocities_medium_complete_peaks; % in AU/s
            matrix_areas_medium_complete_peaks(1:vector_number_medium_complete_peaks(i),i)=vector_areas_medium_complete_peaks; % in AU*s
            
            vector_numbers_large_complete_peaks=find(vector_left_amplitudes_complete_peaks>=(vector_threshold_medium_large_peaks(i)*vector_maximum_left_amplitude_complete_peaks(i)));
            vector_number_large_complete_peaks(i)=length(vector_numbers_large_complete_peaks);
            vector_left_durations_large_complete_peaks=vector_left_durations_complete_peaks(vector_numbers_large_complete_peaks); % in s
            vector_right_durations_large_complete_peaks=vector_right_durations_complete_peaks(vector_numbers_large_complete_peaks); % in s
            vector_left_amplitudes_large_complete_peaks=vector_left_amplitudes_complete_peaks(vector_numbers_large_complete_peaks); % in AU
            vector_left_normalized_amplitudes_large_complete_peaks=vector_left_normalized_amplitudes_complete_peaks(vector_numbers_large_complete_peaks);
            vector_right_amplitudes_large_complete_peaks=vector_right_amplitudes_complete_peaks(vector_numbers_large_complete_peaks); % in AU
            vector_right_normalized_amplitudes_large_complete_peaks=vector_right_normalized_amplitudes_complete_peaks(vector_numbers_large_complete_peaks);
            vector_left_velocities_large_complete_peaks=vector_left_velocities_complete_peaks(vector_numbers_large_complete_peaks); % in AU/s
            vector_right_velocities_large_complete_peaks=vector_right_velocities_complete_peaks(vector_numbers_large_complete_peaks); % in AU/s
            vector_areas_large_complete_peaks=vector_effective_areas_complete_peaks(vector_numbers_large_complete_peaks); % in AU*s
            vector_frequency_large_complete_peaks(i)=vector_number_large_complete_peaks(i)/max(vector_time); % in pk/s
            matrix_numbers_large_complete_peaks(1:vector_number_large_complete_peaks(i),i)=vector_numbers_large_complete_peaks;
            matrix_left_durations_large_complete_peaks(1:vector_number_large_complete_peaks(i),i)=vector_left_durations_large_complete_peaks; % in s
            matrix_right_durations_large_complete_peaks(1:vector_number_large_complete_peaks(i),i)=vector_right_durations_large_complete_peaks; % in s
            matrix_left_amplitudes_large_complete_peaks(1:vector_number_large_complete_peaks(i),i)=vector_left_amplitudes_large_complete_peaks; % in AU
            matrix_left_normalized_amplitudes_large_complete_peaks(1:vector_number_large_complete_peaks(i),i)=vector_left_normalized_amplitudes_large_complete_peaks;            
            matrix_right_amplitudes_large_complete_peaks(1:vector_number_large_complete_peaks(i),i)=vector_right_amplitudes_large_complete_peaks; % in AU
            matrix_right_normalized_amplitudes_large_complete_peaks(1:vector_number_large_complete_peaks(i),i)=vector_right_normalized_amplitudes_large_complete_peaks;
            matrix_left_velocities_large_complete_peaks(1:vector_number_large_complete_peaks(i),i)=vector_left_velocities_large_complete_peaks; % in AU/s
            matrix_right_velocities_large_complete_peaks(1:vector_number_large_complete_peaks(i),i)=vector_right_velocities_large_complete_peaks; % in AU/s
            matrix_areas_large_complete_peaks(1:vector_number_large_complete_peaks(i),i)=vector_areas_large_complete_peaks; % in AU*s
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end

s=menu('Do you want to save the results?','Yes','No');
if s==1
    if operating_system==1
        results_foldername='Results/';
        mkdir([results_foldername rough_data_filename]);
        results_pathname=[results_foldername rough_data_filename '/Results.xlsx'];
    else
        results_foldername='Results\';
        mkdir([results_foldername rough_data_filename]);
        results_pathname=[results_foldername rough_data_filename '\Results.xlsx'];        
    end
    
    for x=1:number_cells
        peak_time=matrix_time_complete_peaks(find(matrix_time_complete_peaks(:,x)),x);
        
        writematrix('peak_time',results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','A1');
        writematrix(peak_time,results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','A2');

        frequency_of_small_Ca_events=vector_frequency_small_complete_peaks(:,x);
        writematrix('frequency_of_small_Ca_events',results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','B1');
        writematrix(frequency_of_small_Ca_events,results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','B2');

        frequency_of_medium_Ca_events=vector_frequency_medium_complete_peaks(:,x);
        writematrix('frequency_of_medium_Ca_events',results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','C1');
        writematrix(frequency_of_medium_Ca_events,results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','C2');

        frequency_of_Ca_transients=vector_frequency_large_complete_peaks(:,x);
        writematrix('frequency_of_Ca_transients',results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','D1');
        writematrix(frequency_of_Ca_transients,results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','D2');

        ascending_time=matrix_left_durations_large_complete_peaks(find(matrix_left_durations_large_complete_peaks(:,x)),x);
        writematrix('ascending_time',results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','E1');
        writematrix(ascending_time,results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','E2');

        decay_time=matrix_right_durations_large_complete_peaks(find(matrix_right_durations_large_complete_peaks(:,x)),x);
        writematrix('decay_time',results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','F1');
        writematrix(decay_time,results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','F2');

        ninety_percent_decay_time=0.9*decay_time;
        writematrix('90%_decay_time',results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','G1');
        writematrix(ninety_percent_decay_time,results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','G2');        
        
        seventy_percent_decay_time=0.7*decay_time;
        writematrix('70%_decay_time',results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','H1');
        writematrix(seventy_percent_decay_time,results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','H2'); 

        fifty_percent_decay_time=0.5*decay_time;
        writematrix('50%_decay_time',results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','I1');
        writematrix(fifty_percent_decay_time,results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','I2'); 
        
        thirty_percent_decay_time=0.3*decay_time;
        writematrix('30%_decay_time',results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','J1');
        writematrix(thirty_percent_decay_time,results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','J2');

        twenty_percent_decay_time=0.2*decay_time;
        writematrix('20%_decay_time',results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','K1');
        writematrix(twenty_percent_decay_time,results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','K2');

        absolute_amplitudes_of_Ca_transients=matrix_left_amplitudes_large_complete_peaks(find(matrix_left_amplitudes_large_complete_peaks(:,x)),x);
        writematrix('absolute_amplitudes_of_Ca_transients',results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','L1');
        writematrix(absolute_amplitudes_of_Ca_transients,results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','L2');

        absolute_amplitude_max=vector_maximum_left_amplitude_complete_peaks(:,x);
        writematrix('absolute_amplitude_max',results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','M1');
        writematrix(absolute_amplitude_max,results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','M2');

        normalised_amplitudes_of_Ca_transients=matrix_left_normalized_amplitudes_large_complete_peaks(find(matrix_left_normalized_amplitudes_large_complete_peaks(:,x)),x);
        writematrix('normalised_amplitudes_of_Ca_transients',results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','N1');
        writematrix(normalised_amplitudes_of_Ca_transients,results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','N2');

        normalised_amplitude_max=vector_maximum_left_normalized_amplitude_complete_peaks(:,x);
        writematrix('normalised_amplitude_max',results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','O1');
        writematrix(normalised_amplitude_max,results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','O2'); 

        velocity=matrix_left_velocities_large_complete_peaks(find(matrix_left_velocities_large_complete_peaks(:,x)),x);
        writematrix('velocity',results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','P1');
        writematrix(velocity,results_pathname,'FileType','spreadsheet','Sheet',['Cell ' num2str(x)],'Range','P2');     
    end
end

close all;
clc;