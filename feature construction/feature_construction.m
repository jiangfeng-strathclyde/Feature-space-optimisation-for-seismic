%% feature construction
clear
clc
close all
%%
load('template noise.mat')
load('template quake.mat')
load('template rockfall.mat')
load('template seismic.mat')

norm_template_quake=template_quake/abs(max(template_quake));
norm_template_seismic=template_seismic/abs(max(template_seismic));
norm_template_rockfall=template_rockfall/abs(max(template_rockfall));
norm_template_noise=template_noise/abs(max(template_noise));


[detected_event,~]=xlsread('detection result');
duration_event=detected_event(:,7);
channel=detected_event(:,10);
detected_time_index=detected_event(:,1:6);
whole_detected_time_index=detected_time_index(:,1)*10^10+detected_time_index(:,2)*10^8+detected_time_index(:,3)*10^6+detected_time_index(:,4)*10^4+detected_time_index(:,5)*60+detected_time_index(:,6);
path=['\...'];
d_file = dir(path);
save_file_name=['feature.mat'];


for i=3:length(d_file)
   
    date_name=d_file(i).name(8:end-4);
    year(i-2)=double(string(date_name(1:4)));
    month(i-2)=double(string(date_name(6:7)));
    day(i-2)=double(string(date_name(9:10)));
    hour(i-2)=double(string(date_name(12:13)));
    minate(i-2)=double(string(date_name(15:16)));
    
    second_sign=string(date_name(end-2));
    
    finger = contains(second_sign,'+');
    isContain=any(finger);
    
    if isContain==1
    second(i-2)=double(string(date_name(18:24)))*10^(double(string(date_name(end-1:end))));
    else 
    second(i-2)=double(string(date_name(18:24)))*10^(-double(string(date_name(end-1:end))));
    end
    
end

whole_true_time_index=year*10^10+month*10^8+day*10^6+hour*10^4+minate*60+second;

for event_index=1:length(whole_detected_time_index)

disp(event_index)
tic
    diff_value=zeros(length(whole_true_time_index),1);
    diff_value = whole_true_time_index-whole_detected_time_index(event_index);
    diff_value = diff_value';
    [min_row,min_col]=find(abs(diff_value)==min(abs(diff_value)));
    time_difference=diff_value(min_row);
    load([path,d_file(min_row+2).name])
    selected_date_name = d_file(min_row+2).name(8:end-4);
    
    file_min=double(string(selected_date_name(15:16)));
    file_second_sign=string(selected_date_name(end-2));
    
    file_finger = contains(file_second_sign,'+');
    isContain=any(file_finger);
    
    if isContain==1
    file_second=double(string(selected_date_name(18:24)))*10^(double(string(selected_date_name(end-1:end))));
    else 
    file_second=double(string(selected_date_name(18:24)))*10^(-double(string(selected_date_name(end-1:end))));
    end
    

    second_Length=2;
    true_event_start=37500-time_difference*250;
    segment=segment;
    detected_event_segment_time_index = [fix(true_event_start-second_Length*250):fix(true_event_start+duration_event(event_index)*250+second_Length*250)];
    detected_event_segment = zeros(length(detected_event_segment_time_index),1);
    
    if true_event_start-second_Length*250<0
        detected_event_segment = [zeros(fix(second_Length*250-true_event_start),1);segment(1:fix(true_event_start+duration_event(event_index)*250+second_Length*250),channel(event_index))];
    elseif true_event_start+duration_event(event_index)*250+second_Length*250>75000
        detected_event_segment = [segment(fix(true_event_start-second_Length*250):75000,channel(event_index));zeros(fix((true_event_start+duration_event(event_index)*250+second_Length*250)-75000),1)];
    else    
    detected_event_segment = segment(detected_event_segment_time_index,channel(event_index));
    end
    
    norm_detected_event_segment = detected_event_segment/(max(detected_event_segment)-min(detected_event_segment));
    
% plor feauture x,y,z
    if true_event_start-second_Length*250<0
        plor_x = [zeros(fix(second_Length*250-true_event_start),1);segment(1:fix(true_event_start+duration_event(event_index)*250+second_Length*250),2)];
        plor_y = [zeros(fix(second_Length*250-true_event_start),1);segment(1:fix(true_event_start+duration_event(event_index)*250+second_Length*250),3)];
        plor_z = [zeros(fix(second_Length*250-true_event_start),1);segment(1:fix(true_event_start+duration_event(event_index)*250+second_Length*250),1)];
    elseif true_event_start+duration_event(event_index)*250+second_Length*250>75000
        plor_x = [segment(fix(true_event_start-second_Length*250):75000,2);zeros(fix((true_event_start+duration_event(event_index)*250+second_Length*250)-75000),1)];
        plor_y = [segment(fix(true_event_start-second_Length*250):75000,3);zeros(fix((true_event_start+duration_event(event_index)*250+second_Length*250)-75000),1)];
        plor_z = [segment(fix(true_event_start-second_Length*250):75000,1);zeros(fix((true_event_start+duration_event(event_index)*250+second_Length*250)-75000),1)];
    else    
    plor_x = segment(detected_event_segment_time_index,2);
    plor_y = segment(detected_event_segment_time_index,3);
    plor_z = segment(detected_event_segment_time_index,1);
    end
    
plor_x=(plor_x-min(plor_x))/(max(plor_x)-min(plor_x));
plor_y=(plor_y-min(plor_y))/(max(plor_y)-min(plor_y));
plor_z=(plor_z-min(plor_z))/(max(plor_z)-min(plor_z));
delt=1;
ttot=length(plor_y);
twin=ttot;

%% signal preprocessing filtering, transfrom 
%setting 
slope_up = []; slope_down=[];slope_up_envelope=[];slope_down_envelope=[];

% filtering
Fs = 250;
fc1 = 5;
fc2 = 10;
fc3 = 50;
fc4 = 70;
fc5 = 100;

[b5_10,a5_10]=butter(4, [fc1 fc2]/(Fs/2));
[b10_50,a10_50]=butter(4, [fc2 fc3]/(Fs/2));
[b5_70,a5_70]=butter(4, [fc1 fc4]/(Fs/2));
[b50_100,a50_100]=butter(4, [fc3 fc5]/(Fs/2));
[b5_100,a5_100]=butter(4, [fc1 fc5]/(Fs/2));

BP_5_10_signal = filter(b5_10,a5_10,norm_detected_event_segment);
BP_10_50_signal = filter(b10_50,a10_50,norm_detected_event_segment);
BP_5_70_signal = filter(b5_70,a5_70,norm_detected_event_segment);
BP_50_100_signal = filter(b50_100,a50_100,norm_detected_event_segment);
BP_5_100_signal = filter(b5_100,a5_100,norm_detected_event_segment);

% envelope signal and calculate the slope
envelope_signal=envelope(norm_detected_event_segment,300); % envelope signal

window_length=20;
for n=2:length(norm_detected_event_segment)/window_length
slope_up(n)=(norm_detected_event_segment(n*window_length)-norm_detected_event_segment(n-1)*window_length)/window_length;
slope_down(n)=(norm_detected_event_segment(n*window_length)-norm_detected_event_segment(n+1)*window_length)/window_length;
slope_up_envelope(n)=(envelope_signal(n*window_length)-envelope_signal(n-1)*window_length)/window_length;
slope_down_envelope(n)=(envelope_signal(n*window_length)-envelope_signal(n+1)*window_length)/window_length;
end

% fourier transform 
N=length(norm_detected_event_segment);
xdft1 = fft(norm_detected_event_segment);
xdft1 = xdft1(1:N/2+1);
xdft1=abs(xdft1); % fft signal
freq = 0:Fs/N:Fs/2;
psdx1 = (1/(N*Fs)) * abs(xdft1).^2; 
psdx1(2:end-1) = 2*psdx1(2:end-1);% power spectral density signal
power_dB=10*log10(psdx1); 

% CWT transfrom
opt.type = 'morlet';         % Mother wavelet type
opt.padtype = 'symmetric';   % padded via symmetrization
opt.rpadded = 1;
opt.nv = 16; 
data.dt=1/Fs;
dt=1/Fs;

data_template_seismic=linspace(0,(dt*length(template_seismic)),length(template_seismic));
[wl_template_seismic,wlas_template_seismic,wldWx_template_seismic] = cwt_fw(template_seismic,opt.type,opt.nv,data.dt);

data_template_noise=linspace(0,(dt*length(template_noise)),length(template_noise));
[wl_template_noise,wlas_template_noise,wldWx_template_noise] = cwt_fw(template_noise,opt.type,opt.nv,data.dt);

data_template_quake=linspace(0,(dt*length(template_quake)),length(template_quake));
[wl_template_quake,wlas_template_quake,wldWx_template_quake] = cwt_fw(template_quake,opt.type,opt.nv,data.dt);

data_template_rockfall=linspace(0,(dt*length(template_rockfall)),length(template_rockfall));
[wl_template_rockfall,wlas_template_rockfall,wldWx_template_rockfall] = cwt_fw(template_rockfall,opt.type,opt.nv,data.dt);

data.x=norm_detected_event_segment;
[wlD,wlasD,wldWxD] = cwt_fw(data.x,opt.type,opt.nv,data.dt);

% discrit wavelet transform
      lev   = 10; wname = 'sym2'; nbcol = 64;
      [signal_c,signal_l] = wavedec(norm_detected_event_segment,lev,wname);
      len = length(norm_detected_event_segment);
      signal_cf = zeros(lev,len);
      d=[];
      for k = 1:lev
         d = detcoef(signal_c,signal_l,k);
         d = d(:)';
         d = d(ones(1,2^k),:);
         signal_cf(k,:) = wkeep1(d(:)',len);
      end
      signal_cf =  signal_cf(:);
      I = find(abs(signal_cf)<sqrt(eps));
      signal_cf(I) = zeros(size(I));
      signal_cf = reshape(signal_cf,lev,len);

     [template_noise_c,template_noise_l] = wavedec(norm_template_noise,lev,wname);
      len = length(template_noise);
      template_noise_cfdD = zeros(lev,len);
      d=[];
      for k = 1:lev
         d = detcoef(template_noise_c,template_noise_l,k);
         d = d(:)';
         d = d(ones(1,2^k),:);
         template_noise_cfdD(k,:) = wkeep1(d(:)',len);
      end
      template_noise_cfdD =  template_noise_cfdD(:);
      I = find(abs(template_noise_cfdD)<sqrt(eps));
      template_noise_cfdD(I) = zeros(size(I));
      template_noise_cfdD = reshape(template_noise_cfdD,lev,len);
      
     [template_quake_c,template_quake_l] = wavedec(norm_template_quake,lev,wname);
      len = length(template_quake);
      template_quake_cfdD = zeros(lev,len);
      d=[];
      for k = 1:lev
         d = detcoef(template_quake_c,template_quake_l,k);
         d = d(:)';
         d = d(ones(1,2^k),:);
         template_quake_cfdD(k,:) = wkeep1(d(:)',len);
      end
      template_quake_cfdD =  template_quake_cfdD(:);
      I = find(abs(template_quake_cfdD)<sqrt(eps));
      template_quake_cfdD(I) = zeros(size(I));
      template_quake_cfdD = reshape(template_quake_cfdD,lev,len);
      
      [template_rockfall_c,template_rockfall_l] = wavedec(norm_template_rockfall,lev,wname);
      len = length(template_rockfall);
      template_rockfall_cfdD = zeros(lev,len);
      d=[];
      for k = 1:lev
         d = detcoef(template_rockfall_c,template_rockfall_l,k);
         d = d(:)';
         d = d(ones(1,2^k),:);
         template_rockfall_cfdD(k,:) = wkeep1(d(:)',len);
      end
      template_rockfall_cfdD =  template_rockfall_cfdD(:);
      I = find(abs(template_rockfall_cfdD)<sqrt(eps));
      template_rockfall_cfdD(I) = zeros(size(I));
      template_rockfall_cfdD = reshape(template_rockfall_cfdD,lev,len);
      
      [template_seismic_c,template_seismic_l] = wavedec(norm_template_seismic,lev,wname);
      len = length(template_seismic);
      template_seismic_cfdD = zeros(lev,len);
      d=[];
      for k = 1:lev
         d = detcoef(template_seismic_c,template_seismic_l,k);
         d = d(:)';
         d = d(ones(1,2^k),:);
         template_seismic_cfdD(k,:) = wkeep1(d(:)',len);
      end
      template_seismic_cfdD =  template_seismic_cfdD(:);
      I = find(abs(template_seismic_cfdD)<sqrt(eps));
      template_seismic_cfdD(I) = zeros(size(I));
      template_seismic_cfdD = reshape(template_seismic_cfdD,lev,len);

  
%%  temporal attributes 
%setting 
duration = []; standerd_deviation=[]; max_envelop_signal=[]; mean_envelop_signal=[]; median_envelop_signal=[]; ratio_mean_max_envelop=[];
ratio_median_max_envelop=[]; rising_duration=[]; decreasing_duration=[]; entropy=[];zcr=[];
skewness_envelop=[]; kurtosis_envelop=[]; rate_of_decay_envelope=[]; ES5_10=[]; ES10_50=[]; ES5_70=[]; ES50_100=[];ES5_100=[];
kurtosis_BP5_10=[]; kurtosis_BP10_50=[]; kurtosis_BP5_70=[]; kurtosis_BP50_100=[]; kurtosis_BP5_100=[]; Energy_average_1_3_autocorr=[]; Energy_average_remaining_autocorr=[];
 num_of_peaks_autocorr=[]; duration_autocorr=[]; RMS_decreasing_l=[]; max_cc_quake=[]; max_cc_rockfall=[]; max_cc_seismic=[]; max_cc_noise=[];
optimum_point_of_separation=[]; ccnAb2_quake=[]; ccnAb2_rockfall=[];ccnAb2_seismic=[];ccnAb2_noise=[];max_cc_noise=[];ypicD_quake=[];ypicD_rockfall=[];ypicD_seismic=[];ypicD_noise=[];
ccnreal2_quake=[]; ccnreal2_rockfall=[]; ccnreal2_seismic=[]; ccnreal2_noise=[]; ccnAb2_quake_without_norm=[]; ccnAb2_rockfall_without_norm=[]; ccnAb2_seismic_without_norm=[];
ccnAb2_noise_without_norm=[]; ccnreal2_quake_without_norm=[]; ccnreal2_rockfall_without_norm=[]; ccnreal2_seismic_without_norm=[]; ccnreal2_noise_without_norm=[];
mean_DFT=[]; max_DFT=[];max_envelop_psd=[]; num_peak_fft=[]; dominant_frequency=[]; spectral_centroid=[]; skewness_power_signal=[]; kurtosis_power_signal=[];
measure_location=[]; measure_of_dispersion=[]; measure_of_assymmetry=[]; measure_of_concentration_around_single_value=[]; first_centbin=[];second_centbin=[];
median_norm_xdft1=[];var_norm_xdft1=[];num_peak_fft_75=[];num_of_peaks_xdft1=[];Energy_0_25_xdft1=[];Energy_25_50_xdft1=[];Energy_50_75_xdft1=[];Energy_75_1_xdft1=[];gyration_radius=[];
sepctral_centroid_width=[];kurtosis_max_stft_over_time=[];mean_ratio_max_mean_stft=[];mean_ratio_max_median_stft=[];num_max_peak_stft_temporal_evolution=[];num_mean_peak_stft_temporal_evolution=[];
num_median_peak_stft_temporal_evolution=[];ratio_82_83=[];ratio_82_84=[];num_peak_stft_central_fre=[];num_peak_stft_max_fre=[];mean_dis_max_mean_fre=[];mean_dis_max_median_fre=[];
mean_distance_first_medain=[];mean_distance_third_medain=[];mean_distance_first_third=[];gamma1=[];gamma2=[];gamma3=[];mean_frequency=[];frequency_bandwidth=[];minimal_frequency=[];
maximal_frequency=[];xpp_quake=[];xpp_rockfall=[];xpp_seismic=[];xpp_noise=[];xpp2D_quake=[];xpp2D_rockfall=[];xpp2D_seismic=[];xpp2D_noise=[];sem_quake_SD=[];sem_rockfall_SD=[];
sem_seismic_SD=[];sem_noise_SD=[];sem_quake_dD=[];sem_seismic_dD=[];sem_rockfall_dD=[];sem_noise_dD=[];std_cepstrum_normal_signal=[];skewness_cepstrum_normal_signal=[];
kurtosis_cepstrum_normal_signal=[];max_cepstrum_normal_signal=[];num_cepstrum_peak=[];num_peaks_lpc=[]; azi=[]; inci=[]; maxeig=[]; Dlp=[]; Dpp=[];

%1. duration
durat = length(detected_event_segment_time_index);
duration = [duration durat];
%2. standerd_deviation
standerd_deviat = std(norm_detected_event_segment)/duration; 
standerd_deviation = [standerd_deviat standerd_deviation];
%3. maxenvelop_signal
max_envelop_sig = max(envelope_signal);
max_envelop_signal = [max_envelop_signal max_envelop_sig];
%4. mean_envelope_signal
mean_envelop_sig = mean(envelope_signal);
mean_envelop_signal = [mean_envelop_signal mean_envelop_sig];
%5. median_envelop_signal
median_envelop_sig = median(envelope_signal);
median_envelop_signal=[median_envelop_signal median_envelop_sig];
%6. ratio of the mean over the maximum of the envelop signal 
ratio_mean_max_enve = mean_envelop_signal/max_envelop_signal;
ratio_mean_max_envelop = [ratio_mean_max_envelop ratio_mean_max_enve];
%7. ratio of the median over the maximum of the envelop signal 
ratio_median_max_enve = median_envelop_signal/max_envelop_signal;
ratio_median_max_envelop = [ratio_median_max_envelop ratio_median_max_enve];
%8. rising_duration
rising_durat=find(norm_detected_event_segment==max(norm_detected_event_segment), 1 );
rising_duration = [rising_duration rising_durat];
%9. decreasing_duration
decreasing_durat=duration-rising_duration;
decreasing_duration = [decreasing_duration decreasing_durat];
%10. entropy norm
entro=yyshang(norm_detected_event_segment,10);
entropy=[entropy entro];
%11. zero_crossing_rate
zc = sum(abs(diff(norm_detected_event_segment>0)))/length(norm_detected_event_segment);  % 11. zeros crossing rate
zcr=[zcr zc];
%12. skewness_envelope
skewness_enve = skewness(envelope_signal);
skewness_envelop =[skewness_envelop skewness_enve];
%13. kurtosis_envelope
kurtosis_enve = kurtosis(envelope_signal);
kurtosis_envelop=[kurtosis_envelop kurtosis_enve];
%14. rate_of_decay_envelope
rate_of_decay_envel=min(slope_down_envelope); 
rate_of_decay_envelope=[rate_of_decay_envelope rate_of_decay_envel];
%15. Energy of the signal filtered in 5–10 Hz
E5_10=trapz([1:N]',envelope(BP_5_10_signal)');
ES5_10 = [ES5_10 E5_10];
%16. Energy of the signal filtered in 10–50 Hz
E10_50=trapz([1:N]',envelope(BP_10_50_signal)');
ES10_50 = [ES10_50 E10_50];
%17. Energy of the signal filtered in 5–70 Hz
E5_70=trapz([1:N]',envelope(BP_5_70_signal)');
ES5_70 = [ES5_70 E5_70];
%18. Energy of the signal filtered in 50–100 Hz
E50_100=trapz([1:N]',envelope(BP_50_100_signal)');
ES50_100 = [ES50_100 E50_100];
%19. Energy of the signal filtered in 5–100 Hz
E5_100=trapz([1:N]',envelope(BP_5_100_signal)');
ES5_100=[ES5_100 E5_100];
%20. Kurtosis of the signal filtered in 5–10 Hz
kurtosis_5_10 = kurtosis(BP_5_10_signal);
kurtosis_BP5_10 = [kurtosis_BP5_10 kurtosis_5_10];
%21. Kurtosis of the signal filtered in 10–50 Hz
kurtosis_10_50 = kurtosis(BP_10_50_signal);
kurtosis_BP10_50=[kurtosis_BP10_50 kurtosis_10_50];
%22. Kurtosis of the signal filtered in 5–70 Hz
kurtosis_5_70 = kurtosis(BP_5_70_signal);
kurtosis_BP5_70=[kurtosis_BP5_70 kurtosis_5_70];
%23. Kurtosis of the signal filtered in 50–100 Hz
kurtosis_50_100 = kurtosis(BP_50_100_signal);
kurtosis_BP50_100=[kurtosis_BP50_100 kurtosis_50_100];
%24. Kurtosis of the signal filtered in 5–100 Hz
kurtosis_5_100 = kurtosis(BP_5_100_signal);
kurtosis_BP5_100=[kurtosis_BP5_100 kurtosis_5_100];
%25. Energy in the first third part of the autocorrelation function
autocorr=[];
autocorr=auto_test(norm_detected_event_segment');
Energy_average_1_3_autoco = sum(abs(autocorr(1:fix(N/3))).^2)/fix(N/3);
Energy_average_1_3_autocorr =[Energy_average_1_3_autocorr Energy_average_1_3_autoco];
%26. Energy in the remaining part of the autocorrelation function
Energy_average_remaining_autoco = sum(abs(autocorr(fix(N/3):end)).^2)/length(autocorr(fix(N/3):end));
Energy_average_remaining_autocorr = [Energy_average_remaining_autocorr Energy_average_remaining_autoco];
%27. Number of peaks in the autocorrelation function
[pks,locs] = findpeaks(autocorr);
line_hi = (1.96)*(1/sqrt(N))+.05;
line_lo = -(1.96)*(1/sqrt(N))-.05;
[wwww,eeee]=find(pks>line_hi);
[num_of_peaks,~]=size(wwww);        % 41. number of peaks of autocorr
num_of_peaks_autocorr=[num_of_peaks_autocorr num_of_peaks];
%28. duration_small0.2__autocorr
row_index=[];
[row_index,~]=find(autocorr<0.2*max(autocorr));
duration_autoco=length(row_index)/N;
duration_autocorr=[duration_autocorr duration_autoco];
%29. RMS between the decreasing part of the signal and  l(t)= Ymax -(ymax/(tf-tmax))t
t=[];
l=[];
Ymax=max(norm_detected_event_segment);
tmax=find(norm_detected_event_segment==max(norm_detected_event_segment));
tf=length(norm_detected_event_segment);
t=1:length(norm_detected_event_segment(rising_duration:end));
l = Ymax -(Ymax/(tf-tmax))*t;
RMS_decreas_l = rms(norm_detected_event_segment(rising_duration:end)-l');
RMS_decreasing_l=[RMS_decreasing_l RMS_decreas_l];
%30. max value crossCorrquake
cross_cor=xcorr(norm_template_quake,norm_detected_event_segment,'none');
[ccor_max, indx]=max(cross_cor);
max_cc_quake=[max_cc_quake ccor_max];
%31. max value crossCorrrockfall
cross_cor=xcorr(norm_template_rockfall,norm_detected_event_segment,'none');
[ccor_max, indx]=max(cross_cor);
max_cc_rockfall=[max_cc_rockfall ccor_max];
%32. max value crossCorrseismic
cross_cor=xcorr(norm_template_seismic,norm_detected_event_segment,'none');
[ccor_max, indx]=max(cross_cor);
max_cc_seismic=[max_cc_seismic ccor_max];
%33. max value crossCorrnoise
cross_cor=xcorr(norm_template_noise,norm_detected_event_segment,'none');
[ccor_max, indx]=max(cross_cor);
max_cc_noise=[max_cc_noise ccor_max];
%34. optimum point of separation OTSU
[na n] = size(norm_detected_event_segment); col = 0;
for i = 1:n
col = col + abs(norm_detected_event_segment(:,i));
end
xotsu = otsu(na,col); % The Otsu method for finding the optimum point of separation 
optimum_point_of_separation = [optimum_point_of_separation xotsu];
% add zeros for normxcorr2
if length(norm_detected_event_segment)<7001
    corr_norm_detected_event_segment=[norm_detected_event_segment;zeros(7001-length(norm_detected_event_segment),1)];
end
%35. 2dNormCrossCorAbquake
cD_ab_quake = normxcorr2(abs(norm_template_quake),abs(corr_norm_detected_event_segment));   
xD_ab_quake= mean(mean(abs(cD_ab_quake)));    
ccnAb2_quake=[ccnAb2_quake xD_ab_quake];
%36. 2dNormCrossCorAbRockfall
cD_ab_rockfall = normxcorr2(abs(norm_template_rockfall),abs(corr_norm_detected_event_segment));   
xD_ab_rockfall= mean(mean(abs(cD_ab_rockfall)));    
ccnAb2_rockfall=[ccnAb2_rockfall xD_ab_rockfall];
%37. 2dNormCrossCorAbseismic
cD_ab_seismic = normxcorr2(abs(norm_template_seismic),abs(corr_norm_detected_event_segment));   
xD_ab_seismic= mean(mean(abs(cD_ab_seismic)));    
ccnAb2_seismic=[ccnAb2_seismic xD_ab_seismic];
%38. 2dNormCrossCorAbNoise
cD_ab_noise = normxcorr2(abs(norm_template_noise),abs(corr_norm_detected_event_segment));   
xD_ab_noise= mean(mean(abs(cD_ab_noise)));    
ccnAb2_noise=[ccnAb2_noise xD_ab_noise];
%39. 2dNormCrossCorAbquakepick
[ypeakD_quake, xpeakD_quake] = find(cD_ab_quake==max(cD_ab_quake(:)));
ypicD_quake=[ypicD_quake ypeakD_quake];
%40. 2dNormCrossCorAbRockfallpick
[ypeakD_rockfall, xpeakD_rockfall] = find(cD_ab_rockfall==max(cD_ab_rockfall(:)));
ypicD_rockfall=[ypicD_rockfall ypeakD_rockfall];
%41. 2dNormCrossCorAbseismicpick
[ypeakD_seismic, xpeakD_seismic] = find(cD_ab_seismic==max(cD_ab_seismic(:)));
ypicD_seismic=[ypicD_seismic ypeakD_seismic];
%42. 2dNormCrossCorAbNoisepick
[ypeakD_noise, xpeakD_noise] = find(cD_ab_noise==max(cD_ab_noise(:)));
ypicD_noise=[ypicD_noise ypeakD_noise];
%43. 2dNormCrossCorrealquake
cD_real_quake = normxcorr2(real(norm_template_quake),real(corr_norm_detected_event_segment));   
xD_real_quake= mean(mean(abs(cD_real_quake)));    
ccnreal2_quake=[ccnreal2_quake xD_real_quake];
%44. 2dNormCrossCorrealRockfall
cD_real_rockfall = normxcorr2(real(norm_template_rockfall),real(corr_norm_detected_event_segment));   
xD_real_rockfall= mean(mean(abs(cD_real_rockfall)));    
ccnreal2_rockfall=[ccnreal2_rockfall xD_real_rockfall];
%45. 2dNormCrossCorrealseismic
cD_real_seismic = normxcorr2(real(norm_template_seismic),real(corr_norm_detected_event_segment));   
xD_real_seismic= mean(mean(abs(cD_real_seismic)));    
ccnreal2_seismic=[ccnreal2_seismic xD_real_seismic];
%46. 2dNormCrossCorrealNoise
cD_real_noise = normxcorr2(real(norm_template_noise),real(corr_norm_detected_event_segment));   
xD_real_noise= mean(mean(abs(cD_real_noise)));    
ccnreal2_noise=[ccnreal2_noise xD_real_noise];
%47. 2dCrossCorAbquake
norm_template_quake_resampled = interp1(linspace(0,1,length(norm_template_quake)), norm_template_quake, (linspace(0,1,length(norm_detected_event_segment))));
cD_ab_quake_without_norm = corr2(abs(norm_template_quake_resampled),abs(norm_detected_event_segment'));   
xD_ab_quake_without_norm= mean(mean(abs(cD_ab_quake_without_norm)));    
ccnAb2_quake_without_norm =[ccnAb2_quake_without_norm xD_ab_quake_without_norm];
%48. 2dCrossCorAbRockfall
norm_template_rockfall_resampled = interp1(linspace(0,1,length(norm_template_rockfall)), norm_template_rockfall, (linspace(0,1,length(norm_detected_event_segment))));
cD_ab_rockfall_without_norm = corr2(abs(norm_template_rockfall_resampled),abs(norm_detected_event_segment'));   
xD_ab_rockfall_without_norm= mean(mean(abs(cD_ab_rockfall_without_norm)));    
ccnAb2_rockfall_without_norm =[ccnAb2_rockfall_without_norm xD_ab_rockfall_without_norm];
%49. 2dCrossCorAbseismic
norm_template_seismic_resampled = interp1(linspace(0,1,length(norm_template_seismic)), norm_template_seismic, (linspace(0,1,length(norm_detected_event_segment))));
cD_ab_seismic_without_norm = corr2(abs(norm_template_seismic_resampled),abs(norm_detected_event_segment'));   
xD_ab_seismic_without_norm= mean(mean(abs(cD_ab_seismic_without_norm)));    
ccnAb2_seismic_without_norm =[ccnAb2_seismic_without_norm xD_ab_seismic_without_norm];
%50. 2dCrossCorAbNoise
norm_template_noise_resampled = interp1(linspace(0,1,length(norm_template_noise)), norm_template_noise, (linspace(0,1,length(norm_detected_event_segment))));
cD_ab_noise_without_norm = corr2(abs(norm_template_noise_resampled),abs(norm_detected_event_segment'));   
xD_ab_noise_without_norm= mean(mean(abs(cD_ab_noise_without_norm)));    
ccnAb2_noise_without_norm =[ccnAb2_noise_without_norm xD_ab_noise_without_norm];
%51. 2dCrossCorRealquake
cD_real_quake_without_norm = corr2(real(norm_template_quake_resampled),real(norm_detected_event_segment'));   
xD_real_quake_without_norm= mean(mean(abs(cD_real_quake_without_norm)));    
ccnreal2_quake_without_norm =[ccnreal2_quake_without_norm xD_real_quake_without_norm];
%52. 2dCrossCorRealRockfall
cD_real_rockfall_without_norm = corr2(real(norm_template_rockfall_resampled),real(norm_detected_event_segment'));   
xD_real_rockfall_without_norm= mean(mean(abs(cD_real_rockfall_without_norm)));    
ccnreal2_rockfall_without_norm =[ccnreal2_rockfall_without_norm xD_real_rockfall_without_norm];
%53. 2dCrossCorRealseismic
cD_real_seismic_without_norm = corr2(real(norm_template_seismic_resampled),real(norm_detected_event_segment'));   
xD_real_seismic_without_norm= mean(mean(abs(cD_real_seismic_without_norm)));    
ccnreal2_seismic_without_norm =[ccnreal2_seismic_without_norm xD_real_seismic_without_norm];
%54. 2dCrossCorRealNoise
cD_real_noise_without_norm = corr2(real(norm_template_noise_resampled),real(norm_detected_event_segment'));   
xD_real_noise_without_norm= mean(mean(abs(cD_real_noise_without_norm)));    
ccnreal2_noise_without_norm =[ccnreal2_noise_without_norm xD_real_noise_without_norm];
%% Spectral attributes

%55. mean_DFT
me_DFT=mean(abs(xdft1));
mean_DFT = [mean_DFT me_DFT];
%56. max_DFT
m_DFT=max(abs(xdft1));
max_DFT = [max_DFT m_DFT];
%57. max_envelop_psd
envelop_psd=envelope(psdx1);
max_envelop_p=max(envelop_psd);
max_envelop_psd = [max_envelop_psd max_envelop_p];
%58. num_peak_fft
num_peak_FFT_xdft1=[];
num_peak_FFT_xdft1=xdft1;
num_peak_ff=numel(findpeaks(num_peak_FFT_xdft1));
num_peak_fft = [num_peak_fft num_peak_ff];
%59. dominant_frequency
locs=[];
[locs] = find(xdft1==max(xdft1));
dominant_frequen=freq(locs);
dominant_frequency = [dominant_frequency dominant_frequen];
%60. spectral_centroid
spectral_centro=xdft1'*freq'/sum(xdft1);
spectral_centroid=[spectral_centroid spectral_centro];
%61. skewness_power_signal
skewness_power_sign = skewness(power_dB); 
skewness_power_signal = [skewness_power_signal skewness_power_sign];
%62. kurtosis_power_signal
kurtosis_power_sign = kurtosis(power_dB); 
kurtosis_power_signal = [kurtosis_power_signal kurtosis_power_sign];
%63. measure_location
pow = xdft1.*conj(xdft1);
k=1:length(pow);
measure_locati=sum(k*pow);
measure_location = [measure_location measure_locati];
%64. measure_of_dispersion
measure_of_dispers=sqrt(((k-measure_location).^2*pow));
measure_of_dispersion=[measure_of_dispersion measure_of_dispers];
%65. measure_of_assymmetry
measure_of_assymme=((k-measure_location).^3*pow)/(measure_of_dispersion.^3);
measure_of_assymmetry=[measure_of_assymmetry measure_of_assymme];
%66. measure_of_concentration_around_single_value
measure_of_concentration_around_single_va=((k-measure_location).^4*pow)/(measure_of_dispersion.^4);
measure_of_concentration_around_single_value=[measure_of_concentration_around_single_value measure_of_concentration_around_single_va];
%67. Central frequency of the 1st quartile
first_quartile_fft = xdft1(find(xdft1<median(xdft1)));
first_abovecutoff = first_quartile_fft > max(first_quartile_fft) / 2;   %3 dB is factor of 2
first_lowbin  = find(first_abovecutoff, 1, 'first');
first_highbin = sum(first_abovecutoff);
first_cent = sqrt(first_lowbin * first_highbin);   %geometric mean
first_centbin=[first_centbin first_cent];
%68. Central frequency of the 2nd quartile
second_quartile_fft = xdft1(find(median(xdft1)<xdft1<mean(xdft1)));
second_abovecutoff = second_quartile_fft > max(second_quartile_fft) / 2;   %3 dB is factor of 2
second_lowbin  = find(second_abovecutoff, 1, 'first');
second_highbin = sum(second_abovecutoff);
second_cent = sqrt(second_lowbin * second_highbin);   %geometric mean
second_centbin=[second_centbin second_cent];
%69. Median of the normalized DFT
norm_xdft1 = xdft1/max(abs(xdft1));
median_norm_xd = median(norm_xdft1);
median_norm_xdft1=[median_norm_xdft1 median_norm_xd];
%70. Variance of the normalized DFT
var_norm_xd = std(norm_xdft1).^2;
var_norm_xdft1 = [var_norm_xdft1 var_norm_xd];
%71. Number of peaks (>0.75 DFTmax)
a=[];b=[];num_peak_FFT_xdft1=[];
[a,b]=find(xdft1>0.75*max(xdft1));
num_peak_FFT_xdft1=xdft1;
num_peak_FFT_xdft1(b)=0;
num_peak_f75=numel(findpeaks(num_peak_FFT_xdft1));
num_peak_fft_75 = [num_peak_fft_75 num_peak_f75];
%72. Number of peaks in the autocorrelation function
autocorr_xdft1=auto_test(xdft1);
[pks_xdft1,locs_xdft1] = findpeaks(autocorr_xdft1);
N_xdft1 = length(autocorr_xdft1);
line_hi_xdft1 = (1.96)*(1/sqrt(N_xdft1))+.05;
line_lo_xdft1 = -(1.96)*(1/sqrt(N))-.05;
[wwww_xdft1,eeee_xdft1]=find(pks_xdft1>line_hi_xdft1);
[num_of_peaks_xd,~]=size(wwww_xdft1);
num_of_peaks_xdft1 = [num_of_peaks_xdft1 num_of_peaks_xd];
%73. Energy in [0, 1/4] Nyf
Energy_0_25_xd = sum(abs(xdft1(1:fix(N_xdft1/4))).^2);
Energy_0_25_xdft1=[Energy_0_25_xdft1 Energy_0_25_xd];
%74. Energy in [1/4, 1/2] Nyf
Energy_25_50_xd = sum(abs(xdft1(fix(N_xdft1/4):fix(N_xdft1/2))).^2);
Energy_25_50_xdft1=[Energy_25_50_xdft1 Energy_25_50_xd];
%75. Energy in [1/2, 3/4] Nyf
Energy_50_75_xd = sum(abs(xdft1(fix(N_xdft1/2):fix(N_xdft1*3/4))).^2);
Energy_50_75_xdft1=[Energy_50_75_xdft1 Energy_50_75_xd];
%76. Energy in [3/4, 1] Nyf
Energy_75_1_xd = sum(abs(xdft1(fix(N_xdft1*3/4):fix(N_xdft1))).^2);
Energy_75_1_xdft1 = [Energy_75_1_xdft1 Energy_75_1_xd];
%77. Gyration radius
m2 = moment(psdx1,2);
m3 = moment(psdx1,3);
gyration_rad=sqrt(m3/m2);
gyration_radius=[gyration_radius gyration_rad];
%78. Spectral centroid width
sepctral_centroid_wi=sqrt(spectral_centroid^2-gyration_radius^2);
sepctral_centroid_width=[sepctral_centroid_width sepctral_centroid_wi];
%79. Kurtosis of the maximum of all discrete Fourier transforms (DFTs)
nwin = 32;
wind = kaiser(nwin,17);
nlap = nwin-10;
nfft = 256;
a_s = []; a_w = []; a_t = [];
[a_s,a_w,a_t]=spectrogram(norm_detected_event_segment,wind,nlap,nfft,Fs,'yaxis');
max_stft_over_time=max(abs(a_s));
kurtosis_max_stft_over_t=kurtosis(max_stft_over_time);
kurtosis_max_stft_over_time = [kurtosis_max_stft_over_time kurtosis_max_stft_over_t];

%80. mean ratio between the maximum and the mean of all STFTs
mean_stft_over_time=mean(abs(a_s));
mean_ratio_max_mean_s = mean(max_stft_over_time/mean_stft_over_time);
mean_ratio_max_mean_stft = [mean_ratio_max_mean_stft mean_ratio_max_mean_s];

%81. Mean ratio between the maximum and the median of all STFTs
median_stft_over_time=median(abs(a_s));
mean_ratio_max_median_s = median(max_stft_over_time/median_stft_over_time);
mean_ratio_max_median_stft = [mean_ratio_max_median_stft mean_ratio_max_median_s];

%82. Number of peaks in the curve showing the temporal evolution of the DFTs maximum
pks_max=[];locs_max=[];
[pks_max,locs_max] = findpeaks(max_stft_over_time);
num_max_peak_stft_temporal_evolut = length(pks_max);
num_max_peak_stft_temporal_evolution = [num_max_peak_stft_temporal_evolution num_max_peak_stft_temporal_evolut];

%83. Number of peaks in the curve showing the temporal evolution of the DFTs mean
pks_mean=[];locs_mean=[];
[pks_mean,locs_mean] = findpeaks(mean_stft_over_time);
num_mean_peak_stft_temporal_evolut = length(pks_mean);
num_mean_peak_stft_temporal_evolution=[num_mean_peak_stft_temporal_evolution num_mean_peak_stft_temporal_evolut];

%84. Number of peaks in the curve showing the temporal evolution of the DFTs median
pks_median=[];locs_median=[];
[pks_median,locs_median] = findpeaks(median_stft_over_time);
num_median_peak_stft_temporal_evolut = length(pks_median);
num_median_peak_stft_temporal_evolution=[num_median_peak_stft_temporal_evolution num_median_peak_stft_temporal_evolut];

%85. Ratio between 88 and 89
ratio_82_8 = num_max_peak_stft_temporal_evolution/num_mean_peak_stft_temporal_evolution;
ratio_82_83=[ratio_82_83 ratio_82_8];
%86. Ratio between 88 and 90
ratio_82_9 = num_max_peak_stft_temporal_evolution/num_median_peak_stft_temporal_evolution;
ratio_82_84 = [ratio_82_84 ratio_82_9];
%87. Number of peaks in the curve of the temporal evolution of the DFTs central frequency
time_index=[];fre_index=[];
[fre_index, time_index] = size(a_s);
central_fre=zeros(time_index,1);
max_fre=[];
central_fre=[];
frequency_mean=[];
frequency_median=[];
a_spectrum=a_s;
for i=1:time_index
    a_spectrum(abs(a_s(:,i))<0.05*max(abs(a_s(:,i))),i)=0;
    frequency_index_row=[];
    frequency_index_col=[];
    [frequency_index_row,frequency_index_col]=find(abs(a_spectrum(:,i))~=0);
if isempty(frequency_index_row)==1
cut_f1=1;cut_f2=120;
    frequency_mean(i) = 60;
    frequency_median(i) = 60;
else
    cut_f1=a_w(min(frequency_index_row));
    cut_f2=a_w(max(frequency_index_row));
    frequency_mean(i) = mean(a_w(frequency_index_row));
    frequency_median(i) = median(a_w(frequency_index_row));
end
    central_fre(i) = (cut_f1(1)+cut_f2(1))/2;
    max_fre(i) = cut_f2;
end
peak_row=[];peak_col=[];
[peak_row,peak_col] = findpeaks(central_fre);
num_peak_stft_central_f = length(peak_row);
num_peak_stft_central_fre=[num_peak_stft_central_fre num_peak_stft_central_f];
%88. Number of peaks in the curve of the temporal evolution of the DFTs maximum frequency
peak_row_max=[];peak_col_max=[];
[peak_row_max,peak_col_max] = findpeaks(max_fre);
num_peak_stft_max_f = length(peak_row_max);
num_peak_stft_max_fre = [num_peak_stft_max_fre num_peak_stft_max_f];

%89. Mean distance between the curves of the temporal evolution of the DFTs maximum frequency and mean frequency
mean_dis_max_mean_f = mean(frequency_mean-max_fre);
mean_dis_max_mean_fre = [mean_dis_max_mean_fre mean_dis_max_mean_f];
%90. Mean distance between the curves of the temporal evolution of the DFTs maximum frequency and median frequency
mean_dis_max_median_f = mean(frequency_median-max_fre);
mean_dis_max_median_fre = [mean_dis_max_median_fre mean_dis_max_median_f];
%91. Mean distance between the 1st quartile and the median of all DFTs as a function of time
first_quartile_index = fix(length(a_t)/4);
second_quartile_index = 2*fix(length(a_t)/4);
third_quartile_index = 3*fix(length(a_t)/4);
time_index_to_the_end = 1:length(a_t);
time_index_median = fix(median(time_index_to_the_end));
mean_distance_first_med = mean(abs(abs(a_s(:,first_quartile_index))-abs(a_s(:,time_index_median))));
mean_distance_first_medain =[mean_distance_first_medain mean_distance_first_med];
%92. Mean distance between the 3rd quartile and the median of all DFTs as a function of time
mean_distance_third_med = mean(abs(abs(a_s(:,third_quartile_index))-abs(a_s(:,time_index_median))));
mean_distance_third_medain = [mean_distance_third_medain mean_distance_third_med];
%93. Mean distance between the 3rd quartile and the 1st quartile of all DFTs as a function of time
mean_distance_first_th = mean(abs(abs(a_s(:,first_quartile_index))-abs(a_s(:,third_quartile_index))));
mean_distance_first_third = [mean_distance_first_third mean_distance_first_th];
%94. gamma1
gamma1_a=freq*(xdft1.^2);
gamma1_b=sum(xdft1.^2);
gamm=gamma1_a/gamma1_b;
gamma1=[gamma1 gamm];
%95. gamma2
gamma2_a=freq.^2*(xdft1.^2);
gamma2_b=sum(xdft1.^2);
gamm=sqrt(gamma2_a/gamma2_b);
gamma2=[gamma2 gamm];
%96. gamma3
gamm=sqrt(abs(gamma1^2-gamma2^2));
gamma3=[gamma3 gamm];
%97. mean frequency
psdx1=psdx1';
mean_frequen=(psdx1*freq')/sum(psdx1);
mean_frequency=[mean_frequency mean_frequen];
%98. frequency_bandwidth
freq_power=freq.^2;
frequency_bandwid=2*sqrt((psdx1*freq_power')/sum(psdx1)-mean_frequency^2);
frequency_bandwidth=[frequency_bandwidth frequency_bandwid];
%99. minimal_frequency
fre_row=[];fre_col=[];
[fre_row,fre_col]=find(psdx1>0.2*max(psdx1));
minimal_frequen=freq(fre_col(1));
minimal_frequency = [minimal_frequency minimal_frequen];
%100. maximal_frequency
maximal_frequen=freq(fre_col(end));
maximal_frequency = [maximal_frequency maximal_frequen];
%101. 2dCrossCorAbquakeDWT
a=[];b=[];
[a,b]=size(signal_cf);
signal_cf_resample=zeros(a,b);
template_quake_cfdD_resample=zeros(a,b);
template_rockfall_cfdD_resample=zeros(a,b);
template_seismic_cfdD_resample=zeros(a,b);
template_noise_cfdD_resample=zeros(a,b);
for i =1:a
template_quake_cfdD_resample(i,:)=interp1(linspace(0,1,length(template_quake_cfdD)), template_quake_cfdD(i,:), (linspace(0,1,length(signal_cf))));
template_rockfall_cfdD_resample(i,:)=interp1(linspace(0,1,length(template_rockfall_cfdD)), template_rockfall_cfdD(i,:), (linspace(0,1,length(signal_cf))));
template_seismic_cfdD_resample(i,:)=interp1(linspace(0,1,length(template_seismic_cfdD)), template_seismic_cfdD(i,:), (linspace(0,1,length(signal_cf))));
template_noise_cfdD_resample(i,:)=interp1(linspace(0,1,length(template_noise_cfdD)), template_noise_cfdD(i,:), (linspace(0,1,length(signal_cf))));
end

cx_quake = normxcorr2(signal_cf,template_quake_cfdD_resample);
xx_quake= mean(mean(abs(cx_quake)));
xpp_quake=[xpp_quake xx_quake];    
   
%102. 2dCrossCorAbRockfallDWT
cx_rockfall = normxcorr2(signal_cf,template_rockfall_cfdD_resample);
xx_rockfall= mean(mean(abs(cx_rockfall)));
xpp_rockfall=[xpp_rockfall xx_rockfall];

%103. 2dCrossCorAbseismicDWT
cx_seismic = normxcorr2(signal_cf,template_seismic_cfdD_resample);
xx_seismic= mean(mean(abs(cx_seismic)));
xpp_seismic=[xpp_seismic xx_seismic];

%104. 2dCrossCorAbNoiseDWT
cx_noise = normxcorr2(signal_cf,template_noise_cfdD_resample);
xx_noise= mean(mean(abs(cx_noise)));
xpp_noise=[xpp_noise xx_noise];

%105. 2dCrossCorRealquakeDWT
cx_quake = corr2(real(signal_cf),real(template_quake_cfdD_resample));
xpp2D_quake=[xpp2D_quake cx_quake];

%106. 2dCrossCorRealRockfallDWT
cx_rockfall = corr2(real(signal_cf),real(template_rockfall_cfdD_resample));
xpp2D_rockfall=[xpp2D_rockfall cx_rockfall];

%107. 2dCrossCorRealseismicDWT
cx_seismic = corr2(real(signal_cf),real(template_seismic_cfdD_resample));
xpp2D_seismic=[xpp2D_seismic cx_seismic];

%108. 2dCrossCorRealNoiseDWT
cx_noise = corr2(real(signal_cf),real(template_noise_cfdD_resample));
xpp2D_noise=[xpp2D_noise cx_noise];



%% cepstrum

%109. std_cepstrum_normal_signal
cepstrum_normal_signal = rceps(norm_detected_event_segment);
length_cepstrum_normal_signal=length(cepstrum_normal_signal);
cepstrum_normal_signal=cepstrum_normal_signal(1:fix(length_cepstrum_normal_signal/2));
std_cepstrum_normal_sig=std(cepstrum_normal_signal);
std_cepstrum_normal_signal=[std_cepstrum_normal_signal std_cepstrum_normal_sig];

%110. skewness_cepstrum_normal_signal
skewness_cepstrum_normal_sig=skewness(cepstrum_normal_signal);
skewness_cepstrum_normal_signal=[skewness_cepstrum_normal_signal skewness_cepstrum_normal_sig];
%111. kurtosis_cepstrum_normal_signal
kurtosis_cepstrum_normal_sig=kurtosis(cepstrum_normal_signal);
kurtosis_cepstrum_normal_signal=[kurtosis_cepstrum_normal_signal kurtosis_cepstrum_normal_sig];
%112. max_cepstrum_normal_signal
max_cepstrum_normal_sig=max(cepstrum_normal_signal);
max_cepstrum_normal_signal=[max_cepstrum_normal_signal max_cepstrum_normal_sig];
%113. cepstrum number of picks (echo)
px=[];locs=[];
[px,locs] = findpeaks(cepstrum_normal_signal,'Threshold',0.2,'MinPeakDistance',0.2);
num_cepstrum_p = length(px);
num_cepstrum_peak=[num_cepstrum_peak num_cepstrum_p];
%% LPC
%114. num of peaks LPC 
lpc_a=[];
lpc_a = lpc(norm_detected_event_segment,8);
signal_est_lpc = filter([0 -lpc_a(2:end)],1,norm_detected_event_segment);
[px_lpc,locs_lpc] = findpeaks(signal_est_lpc,'Threshold',0.1,'MinPeakDistance',0.1);
num_peaks_l = length(px_lpc);
num_peaks_lpc=[num_peaks_lpc num_peaks_l];
% polarity features
[azi,inci,maxeig,Dlp,Dpp] = polarize_estimation(plor_x,plor_y,plor_z,delt,ttot,twin);
toc
event_year=[];event_month=[];event_day=[];event_hour=[];event_min=[];event_sec=[];
event_year=[event_year detected_time_index(event_index,1)];
event_month=[event_month detected_time_index(event_index,2)];
event_day=[event_day detected_time_index(event_index,3)];
event_hour=[event_hour detected_time_index(event_index,4)];
event_min=[event_min detected_time_index(event_index,5)];
event_sec=[event_sec detected_time_index(event_index,6)];

constructed_feature(event_index,1:132)=[event_year event_month event_day event_hour event_min event_sec duration standerd_deviation max_envelop_signal ...  
mean_envelop_signal  median_envelop_signal  ratio_mean_max_envelop ratio_median_max_envelop  rising_duration ...
decreasing_duration entropy zcr skewness_envelop ...
kurtosis_envelop rate_of_decay_envelope  ES5_10  ES10_50 ES5_70 ES50_100 ES5_100 ...
kurtosis_BP5_10  kurtosis_BP10_50  kurtosis_BP5_70  kurtosis_BP50_100  kurtosis_BP5_100  Energy_average_1_3_autocorr  Energy_average_remaining_autocorr ...  
 num_of_peaks_autocorr duration_autocorr  RMS_decreasing_l  max_cc_quake  max_cc_rockfall  max_cc_seismic  max_cc_noise ...
optimum_point_of_separation  ccnAb2_quake  ccnAb2_rockfall ccnAb2_seismic ccnAb2_noise ypicD_quake ypicD_rockfall ypicD_seismic ypicD_noise ...
ccnreal2_quake  ccnreal2_rockfall  ccnreal2_seismic  ccnreal2_noise  ccnAb2_quake_without_norm  ccnAb2_rockfall_without_norm  ccnAb2_seismic_without_norm ... 
ccnAb2_noise_without_norm  ccnreal2_quake_without_norm  ccnreal2_rockfall_without_norm  ccnreal2_seismic_without_norm  ccnreal2_noise_without_norm ...
mean_DFT  max_DFT max_envelop_psd  num_peak_fft  dominant_frequency  spectral_centroid  skewness_power_signal  kurtosis_power_signal ...
measure_location  measure_of_dispersion  measure_of_assymmetry  measure_of_concentration_around_single_value first_centbin second_centbin ... 
median_norm_xdft1 var_norm_xdft1 num_peak_fft_75 num_of_peaks_xdft1 Energy_0_25_xdft1 Energy_25_50_xdft1 Energy_50_75_xdft1 Energy_75_1_xdft1 gyration_radius ... 
sepctral_centroid_width kurtosis_max_stft_over_time mean_ratio_max_mean_stft mean_ratio_max_median_stft num_max_peak_stft_temporal_evolution num_mean_peak_stft_temporal_evolution ...
num_median_peak_stft_temporal_evolution ratio_82_83 ratio_82_84 num_peak_stft_central_fre num_peak_stft_max_fre mean_dis_max_mean_fre mean_dis_max_median_fre ...
mean_distance_first_medain mean_distance_third_medain mean_distance_first_third gamma1 gamma2 gamma3 mean_frequency frequency_bandwidth minimal_frequency ...
maximal_frequency xpp_quake xpp_rockfall xpp_seismic xpp_noise xpp2D_quake xpp2D_rockfall xpp2D_seismic xpp2D_noise  ...
std_cepstrum_normal_signal skewness_cepstrum_normal_signal ...
kurtosis_cepstrum_normal_signal max_cepstrum_normal_signal num_cepstrum_peak num_peaks_lpc azi inci maxeig Dlp Dpp]; 

end


