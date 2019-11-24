% AF_detection_EMD.m

% Copyright (c) Prasanth "Prash" Ganesan
% Author email: <prasganesan.pg@gmail.com>

clear; clc; close all

%% Import AF ECG
[t_af,ECG_af]=rdsamp('afdb/04043',[1],1000); %04015
fs=250;
time_req=3; % required time is 3 sec
samples = time_req*fs;
new_time = linspace(0,time_req,samples);
new_ecg_AF=ECG_af(1:samples);
figure; plot(new_time,new_ecg_AF);
title('3 sec of Atrial Fibrillation Original ECG')
xlabel('Time (sec)')
ylabel('Amplitude')

%% FFT of AF ECG
% Plot in Frequency domain
N_fft= length(new_ecg_AF);
L = length(new_ecg_AF);
new_ecg_AF_fft = abs(fft(new_ecg_AF));
f = [0:(N_fft-1)]*fs/L;
figure; plot(f,abs(new_ecg_AF_fft))
xlabel('frequency (Hz)'); title('Frequqncy domain of AF ECG')

%% Bandpass butterworth filter
[b,a] = butter(1,[0.5 45]./(2*fs),'bandpass');
filtered_new_ecg_AF = filter(b,a,new_ecg_AF);
figure; plot(new_time,filtered_new_ecg_AF);
title('3 sec of Atrial Fibrillation Filtered ECG')
xlabel('Time (sec)')
ylabel('Amplitude')

%% Find annotations
[pkamp,pktime_AF] = findpeaks(filtered_new_ecg_AF,'MinPeakHeight',0.5);
% Plot annotations over that
hold on;
plot(new_time(pktime_AF),pkamp,'r*');

%% Import NSR ecg
[t_nsr,ECG_nsr]=rdsamp('nsrdb/16265',[1],1000);
ann_nsr=rdann('nsrdb/16265','atr',[],1000);
fs=128;
time_req=3; % required time is 3 sec
samples = time_req*fs;
new_time = linspace(0,time_req,samples);
new_ecg_NSR=ECG_nsr(1:samples);
figure; plot(new_time,new_ecg_NSR);
title('3 sec of NSR Original ECG')
xlabel('Time (sec)')
ylabel('Amplitude')

%% FFT of NSR ECG
% Plot in Frequency domain
N_fft= length(new_ecg_NSR);
L = length(new_ecg_NSR);
new_ecg_NSR_fft = abs(fft(new_ecg_NSR));
f = [0:(N_fft-1)]*fs/L;
figure; plot(f,abs(new_ecg_NSR_fft))
xlabel('frequency (Hz)'); title('Frequqncy domain of NSR ECG')

%% Bandpass butterworth filter
[b,a] = butter(1,[0.5 45]./(2*fs),'bandpass');
filtered_new_ecg_NSR = filter(b,a,new_ecg_NSR);
figure; plot(new_time,filtered_new_ecg_NSR);
title('3 sec of NSR Filtered ECG')
xlabel('Time (sec)')
ylabel('Amplitude')

%% Find annotations for NSR
[pkamp,pktime_NSR] = findpeaks(filtered_new_ecg_NSR,'MinPeakHeight',0.5);
% Plot annotations over that
hold on;
plot(new_time(pktime_NSR),pkamp,'r*');

%% Choose One R-R cycle
cycle_AF = filtered_new_ecg_AF(pktime_AF(2):pktime_AF(3));
cycle_NSR = filtered_new_ecg_NSR(pktime_NSR(2):pktime_NSR(3));
close all
figure; plot(cycle_AF);
title('One cycle of AF')
xlabel('Samples')
ylabel('Amplitude')

figure; plot(cycle_NSR);
title('One cycle of NSR')
xlabel('Samples')
ylabel('Amplitude')

%% Apply EMD
close all
[imf_AF,~,~] = emd(cycle_AF);
[imf_NSR,~,~] = emd(cycle_NSR);
% figure;
for i=1:size(imf_NSR,1)
%     subplot(size(imf_NSR,1),1,i)
    figure;
    plot(imf_NSR(i,:))
    title(['NSR IMF ' num2str(i)])
    xlabel('samples')
    ylabel('amplitude')
end

% figure;
for i=1:size(imf_AF,1)
%     subplot(size(imf_AF,1),1,i)
    figure;
    plot(imf_AF(i,:))
    title(['AF IMF ' num2str(i)])
    xlabel('samples')
    ylabel('amplitude')
end

%% Plot the positive halves of the IMF3 and IMF4
% close all
% imf3_AF = imf_AF(3,:);
% imf3_NSR = imf_NSR(3,:);
% imf4_AF = imf_AF(4,:);
% imf4_NSR = imf_NSR(4,:);
% 
% figure;
% plot(imf3_AF(imf3_AF>0))
% title('AF IMF 3 positive')
% 
% figure;
% plot(imf3_NSR(imf3_NSR>0))
% title('NSR IMF 3 positive')
% 
% figure;
% plot(imf4_AF(imf4_AF>0))
% title('AF IMF 4 positive')
% 
% figure;
% plot(imf4_NSR(imf4_NSR>0))
% title('NSR IMF 4 positive')

%% Statistical ANalysis in EMD DOmain
close all
% Variance and mean of AF IMF
for i=1:size(imf_AF,1)
variances_AF_IMF(i) = var(imf_AF(i,:));
end

for i=1:size(imf_AF,1)
means_AF_IMF(i) = mean(imf_AF(i,:));
end

% Variance and mean of NSR IMF
for i=1:size(imf_NSR,1)
variances_NSR_IMF(i) = var(imf_NSR(i,:));
end

for i=1:size(imf_NSR,1)
means_NSR_IMF(i) = mean(imf_NSR(i,:));
end

figure; plot(variances_AF_IMF); title('Variances of AF');
xlabel('IMF No.'); ylabel('Variance')
figure; plot(variances_NSR_IMF); title('Variances of NSR');
xlabel('IMF No.'); ylabel('Variance')
figure; plot(means_AF_IMF); title('Means of AF');
xlabel('IMF No.'); ylabel('Mean')
figure; plot(means_NSR_IMF); title('Means of NSR');
xlabel('IMF No.'); ylabel('Mean')
