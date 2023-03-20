clc;
clear;
close all;

%% Question 2

% Loading the needed variables
load("time.mat","time")
load("voltage1.mat","voltage1")
load("voltage2.mat","voltage2")
load("voltage3.mat","voltage3")
load("current.mat","current")

% FFT needed variables
NFFT = 20000;
L = length(time);
Fs = 10^5;

% Fourier Transform
f= Fs/2*linspace(0,1,NFFT/2+1);

% Calculating the fourier of each signal
Y_vol1 = fft(voltage1,NFFT)/L;
Y_vol2 = fft(voltage2,NFFT)/L;
Y_vol3 = fft(voltage3,NFFT)/L;
Y_curr = fft(current,NFFT)/L;


%% Question 3

% Exctracting half of the signal
Y_vol1_h = 2*abs(Y_vol1(1:NFFT/2+1));
Y_vol2_h = 2*abs(Y_vol2(1:NFFT/2+1));
Y_vol3_h = 2*abs(Y_vol3(1:NFFT/2+1));
Y_curr_h = 2*abs(Y_curr(1:NFFT/2+1));



% Limiting the frequency axis
FREQ = 1000;

% Figure fixed size
size = [10 10 1300 600];


%%%%%%%%%%%%%%%%%%%%%%55
% Plotting current along side its spectral content
f1 = figure('Renderer', 'painters','Name','Voltage_1','NumberTitle','off', 'Position', size);

% Plotting the current signal in time
subplot(1,2,1);
plt_volt1_T = plot(time, voltage1);
xlabel('time (sec)');
ylabel('Voltage_1(V)');
title('Voltage_1 signal in time')
axis([0 ,max(time), -150, 150])
% Data tips
% target_value = 0.1;
% indx = find(time<=target_value, 1, 'last');
% dt_curr_t = datatip(plt_curr_T, target_value, current(indx));
% dt_curr_t.Location = 'southeast';

% Plotting the single sided spectral content of current 
subplot(1,2,2);
plt_volt1_F = plot(f,Y_vol1_h);
xlim([0 FREQ])
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');
title('Single Sided Amplituted Spectrum of Voltage_1')


% Saving the figure
print -depsc voltage1_Q3


% Plotting current along side its spectral content
f4 = figure('Renderer', 'painters','Name','Current','NumberTitle','off', 'Position', [10 10 1300 600]);

% Plotting the current signal in time
subplot(1,2,1);
plt_curr_T = plot(time, current);
xlabel('time (sec)');
ylabel('Current(t)');
title('Current signal in time')
% Data tips
target_value = 0.1;
indx = find(time<=target_value, 1, 'last');
dt_curr_t = datatip(plt_curr_T, target_value, current(indx));
dt_curr_t.Location = 'southeast';

% Plotting the single sided spectral content of current 
subplot(1,2,2);
plt_curr_F = plot(f,Y_curr_h);
xlim([0 FREQ])
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');
title('Single Sided Amplituted Spectrum of Current')
% Data tips
target_value = 1/target_value;
indx = find(f<=target_value, 1, 'last');
dt_curr_t = datatip(plt_curr_F, target_value, Y_curr_h(indx));
dt_curr_t.Location = 'northeast';

% Saving the figure
print -depsc current_Q3







% % Plotting each fourier
% f11 = figure('Renderer', 'painters','Name','FFTs of Signals','NumberTitle','off');
% 
% % Current
% subplot(2,2,1);
% plot(f,Y_curr_h)
% xlabel('Frequency (Hz)');
% ylabel('|Y(f)|');
% title('Fourrier transform Current')
% 
% % Voltage1
% subplot(2,2,2);
% plot(f,Y_vol1_h)
% title('Fourrier transform voltage1');
% 
% % Voltage2
% subplot(2,2,3);
% plot(f,Y_vol2_h)
% title('Fourrier transform voltage2');
% 
% % Voltage3
% subplot(2,2,4);
% plot(f,Y_vol3_h)
% title('Fourrier transform voltage3');
% 
% 
% %% Question 4
% 
% %..................................THD.....................................
% % Extracting the fundamental harmonic index
% nominal_frequency = 10;
% fundamental_harmonic_indx = find(f == nominal_frequency);
% 
% % Current
% fundamental_amplitude_curr = Y_curr_h(fundamental_harmonic_indx);
% curr_rms = fundamental_amplitude_curr/sqrt(2);
% harmonic_pw = sum((Y_curr_h / sqrt(2)).^2) - (curr_rms)^2;
% THD_curr = 100 * sqrt(harmonic_pw) / curr_rms;
% 
% % voltage1
% fundamental_amplitude_vol1 = Y_vol1_h(fundamental_harmonic_indx);
% volt1_rms = fundamental_amplitude_vol1/sqrt(2);
% harmonic_pw = sum((Y_vol1_h / sqrt(2)).^2) - (volt1_rms)^2;
% THD_vol1 = 100 * sqrt(harmonic_pw) / volt1_rms;
% 
% % voltage2
% fundamental_amplitude_vol2 = Y_vol2_h(fundamental_harmonic_indx);
% volt2_rms = fundamental_amplitude_vol2/sqrt(2);
% harmonic_pw = sum((Y_vol2_h / sqrt(2)).^2) - (volt2_rms)^2;
% THD_vol2 = 100 * sqrt(harmonic_pw) / volt2_rms;
% 
% % voltage3
% fundamental_amplitude_vol3 = Y_vol3_h(fundamental_harmonic_indx);
% volt3_rms = fundamental_amplitude_vol3/sqrt(2);
% harmonic_pw = sum((Y_vol3_h / sqrt(2)).^2) - (volt3_rms)^2;
% THD_vol3 = 100 * sqrt(harmonic_pw) / volt3_rms;
% 
% 
% % ..................................DPF.....................................
% % Calculating the angle of the fundamental harmonic of each signal
% angl_curr = angle(Y_curr(fundamental_harmonic_indx));
% angl_vol1 = angle(Y_vol1(fundamental_harmonic_indx));
% angl_vol2 = angle(Y_vol2(fundamental_harmonic_indx));
% angl_vol3 = angle(Y_vol3(fundamental_harmonic_indx));
% 
% % Calculating the DPF values
% DPF_1 = cos(angl_vol1 - angl_curr);
% DPF_2 = cos(angl_vol2 - angl_curr);
% DPF_3 = cos(angl_vol3 - angl_curr);
% 
% 
% % ...................................PF.....................................
% % Calculating the rms values
% I1_rms = (fundamental_amplitude_curr/sqrt(2));
% I_rms = sum(Y_curr_h/sqrt(2));
% 
% % Calculating the PF values
% PF_1 = I1_rms / I_rms * DPF_1;
% PF_2 = I1_rms / I_rms * DPF_2;
% PF_3 = I1_rms / I_rms * DPF_3;
% 
% 
% %..............................Apparent Power..............................
% % Calculating the rms values
% V1_rms = sum(Y_vol1_h/sqrt(2));
% V2_rms = sum(Y_vol2_h/sqrt(2));
% V3_rms = sum(Y_vol3_h/sqrt(2));
% 
% % Calculating the S values
% S1 = V1_rms * I_rms;
% S2 = V2_rms * I_rms;
% S3 = V3_rms * I_rms;
% 
% %...............................Active Power...............................
% % Calculating the P values
% P1 = S1 * PF_1;
% P2 = S2 * PF_2;
% P3 = S3 * PF_3;
% 
% 
% %% Question 5
% % Calculating the V2_rms value through theory
% V2_rms_theoretical = sqrt(1/length(time) * sum(voltage2.^2));
% 
% % Extracting the index of 700Hz
% sample_indx = find(f==700);
% 
% % Extracting the rms value of Voltage2 througth the values until 700Hz
% V2_rms_experimental = sqrt(sum((Y_vol2_h(1:sample_indx) /sqrt(2)).^2));
