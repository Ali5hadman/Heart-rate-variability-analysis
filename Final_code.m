%% Clear workspace
clear;
close all;
clc;

%% Extracting active and quite sleep

i=5;

load DATA;
fs = 500;

b=LPT_ALLDATA.Ecg(i);
b=b{1};

c=LPT_ALLDATA.SS_baseline(i);
f=cell2mat(c{1});

d_1 = find(f==1);
d_2 = find(f==2);

tt_active=f(d_2-2:d_2-1)*fs;
tt_quite=f(d_1-2:d_1-1)*fs;

if i==2
    tt_active=f(1:2)*fs;
    tt_quite=([185 ,545])*fs;
    active= b(1:tt_active(2));
    quite= b(tt_quite(1):tt_quite(2));
elseif i==5    
    active= b(1:tt_active(2));
    quite= b(tt_quite(1):tt_quite(2));
else
    active= b(tt_active(1):tt_active(2));
    quite= b(tt_quite(1):tt_quite(2));
end

t_active = 0:1/fs:length(active)/fs-1/fs;
t_quite = 0:1/fs:length(quite)/fs-1/fs;


% Geting the spectrum:

active = active-mean(active);
quite = quite-mean(quite);

active_spectrum = abs(fft(active,2048));
quite_spectrum = abs(fft(quite,2048));
fa = linspace(0,fs,length(active_spectrum));
w = fa/fs;



% Filtering:

filter_order = 5 ;
fcut = 0.1;       % cut-off frequency [Hz]
wc = fcut/(fs/2); % NORMALIZED cut-off frequency
[b,a] = butter(filter_order,wc,'high');
figure; freqz(b,a,1024,fs);

active_ECG_filtered = filtfilt(b,a,active);
quite_ECG_filtered = filtfilt(b,a,quite);

activefilt = movmean(active_ECG_filtered,filter_order);
quitefilt = movmean(quite_ECG_filtered,filter_order);



% Finding peaks and RR intervals:

if i==2 
    [Rpeak_active ,loc_a ] = findpeaks(activefilt, t_active, "MinPeakHeight", 3*10^3 , "MinPeakDistance",0.25); 
    [ Rpeak_quite, loc_q] = findpeaks(quitefilt, t_quite, "MinPeakHeight", 3*10^3 , "MinPeakDistance",0.25);
elseif i==3
    [Rpeak_active ,loc_a ] = findpeaks(activefilt, t_active, "MinPeakHeight", 6*10^3 , "MinPeakDistance",0.25); 
    [ Rpeak_quite, loc_q] = findpeaks(quitefilt, t_quite, "MinPeakHeight", 6 *10^3 , "MinPeakDistance",0.25);  
elseif i==4
    [Rpeak_active ,loc_a ] = findpeaks(activefilt, t_active, "MinPeakHeight", 1*10^3 , "MinPeakDistance",0.25); 
    [ Rpeak_quite, loc_q] = findpeaks(quitefilt, t_quite, "MinPeakHeight", 1 *10^3 , "MinPeakDistance",0.25);  
elseif i==5
    [Rpeak_active ,loc_a ] = findpeaks(activefilt, t_active, "MinPeakHeight", 0.8*10^3 , "MinPeakDistance",0.25); 
    [ Rpeak_quite, loc_q] = findpeaks(quitefilt, t_quite, "MinPeakHeight", 0.8*10^3 , "MinPeakDistance",0.25); 
else
    [Rpeak_active ,loc_a ] = findpeaks(activefilt, t_active, "MinPeakHeight", 6*10^3 , "MinPeakDistance",0.25); 
    [ Rpeak_quite, loc_q] = findpeaks(quitefilt, t_quite, "MinPeakHeight", 6*10^3 , "MinPeakDistance",0.25); 
end

RR_activeR = ((diff(loc_a)))*1000; %expressed in msec
RR_quiteR = ((diff(loc_q)))*1000;


% Finding HR and HR mean:

HR_active = 60./(RR_activeR/1000);
HR_quite = 60./(RR_quiteR/1000);

HR_active_mean = mean(HR_active);
HR_quite_mean = mean(HR_quite);

% Resampling for PSD:

%Detrend
RR_active = detrend(RR_activeR);
RR_quite = detrend(RR_quiteR);

%Resampling
f_rs = 10;
RR_active_rs =interp1(loc_a(1:end-1), RR_active,(loc_a:1/f_rs:loc_a(end-1)-1/f_rs));
RR_quite_rs = interp1(loc_q(1:end-1), RR_quite , (loc_q:1/f_rs:loc_q(end-1)-1/f_rs));


% Non-Parametric PSD

if i==3
    window = 100;    %60 samples=30 sec
    overlap = 50;   %overlap 50%
    nfft = 1024;
elseif i==4
    window = 25;    %60 samples=30 sec
    overlap = 13;   %overlap 50%
    nfft = 1024;
elseif i==5
    window = 20;    %60 samples=30 sec
    overlap = 10;   %overlap 50%
    nfft = 1024;
else
    window = 20;    %60 samples=30 sec
    overlap = 10;   %overlap 50%
    nfft = 1024;
end 

% Parametric PSD
AR_order = 8;
[PSD_YW_RR_active,f] = pyulear(RR_active_rs,AR_order,nfft,f_rs);
[PSD_YW_RR_quite,f] = pyulear(RR_quite_rs,AR_order,nfft,f_rs);



%Non-Parametric PSD
[PSD_welch_RR_active,f] = pwelch(RR_active_rs,hamming(window),overlap,nfft,f_rs);
[PSD_welch_RR_quite,f] = pwelch(RR_quite_rs,hamming(window),overlap,nfft,f_rs);



% Finding HF and LF:

%Power indices
LF = [0.05 0.2]; %Low Frequency band
HF = [0.5 1.5];  %High Frequency band


f_LF = and(ge(f,LF(1)),le(f,LF(2))); %ge=greater than or equal to; le=less than or equal to
f_HF = and(ge(f,HF(1)),le(f,HF(2)));

LF_active_welch = trapz(PSD_welch_RR_active(f_LF));
HF_active_welch = trapz(PSD_welch_RR_active(f_HF));
LF_active_YW = trapz(PSD_YW_RR_active(f_LF));
HF_active_YW = trapz(PSD_YW_RR_active(f_HF));

LF_quite_welch = trapz(PSD_welch_RR_quite(f_LF));
HF_quite_welch = trapz(PSD_welch_RR_quite(f_HF));
LF_quite_YW = trapz(PSD_YW_RR_quite(f_LF));
HF_quite_YW = trapz(PSD_YW_RR_quite(f_HF));

LF2HF_active_welch = LF_active_welch/HF_active_welch;
LF2HF_active_YW = LF_active_YW/HF_active_YW;

LF2HF_quite_welch = LF_quite_welch/HF_quite_welch;
LF2HF_quite_YW = LF_quite_YW/HF_quite_YW;


%% Ploting playground:

%% Raw Active and Quite ECG signals:
figure;
subplot(2,1,1); plot(t_active,active);xlim([50 55]); title('Active_ECG', 'Interpreter', 'none'); xlabel('Time [s]'); ylabel('Amplitude [mV]');
subplot(2,1,2); plot(t_quite,quite);xlim([50 55]); title('Quite_ECG', 'Interpreter', 'none'); xlabel('Time [s]'); ylabel('Amplitude [mV]');


%% Active and Quite Spectrum 
figure;
subplot(2,1,1); plot(fa,active_spectrum); title('active_ECG_spectrum', 'Interpreter', 'none'); xlabel('Frequency [Hz]'); ylabel('|X(f)|');
subplot(2,1,2); plot(fa,quite_spectrum); title('quiet_ECG_spectrum', 'Interpreter', 'none'); xlabel('Frequency [Hz]'); ylabel('|X(f)|');


%% Raw and Filtered of Active and Quite ECG 
figure;
subplot(2,1,1); plot(t_active,active); xlim([50 55]);hold on; plot(t_active,activefilt); title('active & active_filtered', 'Interpreter', 'none'); xlabel('Time [s]'); ylabel('Amplitude [mV]'); ylim([-10000 12000]);
subplot(2,1,2); plot(t_quite,quite);xlim([0 180]); hold on; plot(t_quite,quitefilt);  title('quite & quite_filtered', 'Interpreter', 'none'); xlabel('Time [s]'); ylabel('Amplitude [mV]'); ylim([-10000 12000]);
%%
% figure;
% subplot(2,1,1); plot(t_active,active);xlim([50 55]); hold on; plot(t_active,active_ECG_filtered); title('active & active_filtered', 'Interpreter', 'none'); xlabel('Time [s]'); ylabel('Amplitude [mV]'); ylim([-10000 12000]);
% subplot(2,1,2); plot(t_quite,quite);xlim([50 55]) ;hold on; plot(t_quite,quite_ECG_filtered);  title('quite & quite_filtered', 'Interpreter', 'none'); xlabel('Time [s]'); ylabel('Amplitude [mV]'); ylim([-10000 12000]);


%% Peaks and RR signal for Active and Quite ECG 
figure;
subplot(2,2,3); plot(t_active,activefilt); hold on; plot(loc_a, Rpeak_active,'ok');xlim([50 70]); title('ECG_active_filtered & R_peaks_ECG_active_filtered', 'Interpreter', 'none'); xlabel('Time [s]'); ylabel('Amplitude [mV]');
subplot(2,2,1); plot(t_quite,quitefilt); hold on; plot(loc_q, Rpeak_quite,'ok'); xlim([50 70]); title('ECG_quite_filtered & R_peaks_ECG_quite_filtered', 'Interpreter', 'none'); xlabel('Time [s]'); ylabel('Amplitude [mV]');
subplot(2,2,2); plot(RR_activeR,'-*'); title('RR_quite_ECG_filtered', 'Interpreter', 'none'); xlabel('Beat Number'); ylabel('Duration [ms]'); ylim([300 700]); xlim([0 230]);
subplot(2,2,4); plot(RR_quiteR,'-*'); title('RR_active_ECG_filtered', 'Interpreter', 'none'); xlabel('Beat Number'); ylabel('Duration [ms]'); ylim([300 700]); xlim([0 230]);


%% HF and LF for Active and Quite ECG
figure;
subplot(2,1,1); plot(loc_a(1:end-1),RR_active);
title('RR_ECG_active','Interpreter','none'); xlabel('Time [s]'); ylabel('Duration [ms]'); ylim([-150 100]);

subplot(2,1,2); plot(loc_q(1:end-1),RR_quite);
title('RR_ECG_quite','Interpreter','none'); xlabel('Time [s]'); ylabel('Duration [ms]'); ylim([-150 100]);


figure;
subplot(2,2,1); plot(f,PSD_welch_RR_active);
title('PSD_Welch RR_ECG_active','Interpreter','none'); xlabel('Frequency [Hz]'); ylabel('PSD Welch Method [ms^2/Hz]');
str = ['LF/HF=',num2str(LF2HF_active_welch)]; text(0.6,0.8*max(PSD_welch_RR_active),str,'HorizontalAlignment','left');

subplot(2,2,2); plot(f,PSD_YW_RR_active);
title('PSD_YW RR_ECG_active','Interpreter','none'); xlabel('Frequency [Hz]'); ylabel('PSD Yule-Walker Method [ms^2/Hz]');
str = ['LF/HF=',num2str(LF2HF_active_YW)]; text(0.6,0.8*max(PSD_YW_RR_active),str,'HorizontalAlignment','left');

subplot(2,2,3); plot(f,PSD_welch_RR_quite);
title('PSD_Welch RR_ECG_quite','Interpreter','none'); xlabel('Frequency [Hz]'); ylabel('PSD Welch Method [ms^2/Hz]');
str = ['LF/HF=',num2str(LF2HF_quite_welch)]; text(0.6,0.8*max(PSD_welch_RR_quite),str,'HorizontalAlignment','left');

subplot(2,2,4); plot(f,PSD_YW_RR_quite);
title('PSD_YW RR_ECG_quite','Interpreter','none'); xlabel('Frequency [Hz]'); ylabel('PSD Yule-Walker Method [ms^2/Hz]');
str = ['LF/HF=',num2str(LF2HF_quite_YW)]; text(0.6,0.8*max(PSD_YW_RR_quite),str,'HorizontalAlignment','left');



