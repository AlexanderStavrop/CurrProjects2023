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
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Voltage_1 in time and in Frequency spectrum
f1 = figure('Renderer', 'painters','Name','Voltage_1','NumberTitle','off', 'Position', [10 1400 2000 500]);

% Plotting the current signal in time
subplot(1,2,1);
plt_volt1_T = plot(time, voltage1);
xlabel('Time (sec)');
ylabel('voltage\_1 (V)');
title('voltage\_1 signal in time')
axis([0, 0.2, -150, 150])

% Plotting the single sided spec tral content of current 
subplot(1,2,2);
plt_volt1_F = plot(f,Y_vol1_h);
axis([0 1800 0 90 ])
xticks(0:200:1800);
xlabel('Frequency (Hz)');
ylabel('|Y\_vol1 (f)|');
title('Single Sided Amplituted Spectrum of Voltage\_1')

% Data tips
target_value = 10;
indx = find(f<=target_value, 1, 'last');
dt_volt1_t = datatip(plt_volt1_F, target_value, Y_vol1_h(indx));
dt_volt1_t.Location = 'southwest';
target_value = 200;
indx = find(f<=target_value, 1, 'last');
dt_volt1_F = datatip(plt_volt1_F, target_value, Y_vol1_h(indx));
dt_volt1_F.Location = 'southeast';

% Saving the figure
print -depsc Q3_voltage1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Voltage_2 in time and in Frequency spectrum
f2 = figure('Renderer', 'painters','Name','Voltage_2','NumberTitle','off', 'Position',[10 740 2000 500]);

% Plotting the current signal in time
subplot(1,2,1);
plt_volt2_T = plot(time, voltage2);
xlabel('Time (sec)');
ylabel('voltage\_2 (V)');
title('Voltage\_2 signal in time')
axis([0 , 0.2, -150, 150])

% Plotting the single sided spec tral content of current 
subplot(1,2,2);
plt_volt2_F = plot(f,Y_vol2_h);
axis([0 1000 0 90 ])
xticks(0:200:1000);
xlabel('Frequency (Hz)');
ylabel('|Y\_vol2 (f)|');
title('Single Sided Amplituted Spectrum of Voltage\_2')

% Data tips
target_value = 10;
indx = find(f<=target_value, 1, 'last');
dt_volt2_F = datatip(plt_volt2_F, target_value, Y_vol2_h(indx));
dt_volt2_F.Location = 'southeast';
target_value = 200;
indx = find(f<=target_value, 1, 'last');
dt_volt2_F = datatip(plt_volt2_F, target_value, Y_vol2_h(indx));
dt_volt2_F.Location = 'northeast';

% Saving the figure
print -depsc Q3_voltage2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Voltage_3 in time and in Frequency spectrum
f3 = figure('Renderer', 'painters','Name','Voltage_3','NumberTitle','off', 'Position', [10 143 2000 500]);

% Plotting the current signal in time
subplot(1,2,1);
plt_volt3_T = plot(time, voltage3);
xlabel('Time (sec)');
ylabel('voltage\_3 (V)');
title('Voltage\_3 signal in time')
axis([0 , 0.2, -110, 110])

% Plotting the single sided spec tral content of current 
subplot(1,2,2);
plt_volt3_F = plot(f,Y_vol3_h);
axis([0 1000 0 90 ])
xticks(0:200:1000);
xlabel('Frequency (Hz)');
ylabel('|Y\_vol3 (f)|');
title('Single Sided Amplituted Spectrum of Voltage\_3')
% Data tips
target_value = 10;
indx = find(f<=target_value, 1, 'last');
dt_volt3_F = datatip(plt_volt3_F, target_value, Y_vol3_h(indx));
dt_volt3_F.Location = 'southeast';

% Saving the figure
print -depsc Q3_voltage3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Current in time and in Frequency spectrum
f4 = figure('Renderer', 'painters','Name','Current','NumberTitle','off', 'Position', [10 10 1300 600]);

% Plotting the current signal in time
subplot(1,2,1);
plt_curr_T = plot(time, current);
xlabel('Time (sec)');
ylabel('Current (t)');
title('Current signal in time')

% Plotting the single sided spectral content of current 
subplot(1,2,2);
plt_curr_F = plot(f,Y_curr_h);
xlim([0 1000])
xlabel('Frequency (Hz)');
ylabel('|Y\_curr  (f)|');
title('Single Sided Amplituted Spectrum of Current')
% Data tips
target_value = 10;
indx = find(f<=target_value, 1, 'last');
dt_curr_t = datatip(plt_curr_F, target_value, Y_curr_h(indx));
dt_curr_t.Location = 'northeast';

% Saving the figure
print -depsc Q3_current
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Question 4

%..................................THD.....................................
% Extracting the fundamental harmonic index
nominal_frequency = 10;
fundamental_harmonic_indx = find(f == nominal_frequency);

% Current
fundamental_amplitude_curr = Y_curr_h(fundamental_harmonic_indx);
curr_rms = fundamental_amplitude_curr/sqrt(2);
harmonic_pw = sum((Y_curr_h / sqrt(2)).^2) - (curr_rms)^2;
THD_curr = 100 * sqrt(harmonic_pw) / curr_rms;

% voltage_1
fundamental_amplitude_vol1 = Y_vol1_h(fundamental_harmonic_indx);
volt1_rms = fundamental_amplitude_vol1/sqrt(2);
harmonic_pw = sum((Y_vol1_h / sqrt(2)).^2) - (volt1_rms)^2;
THD_vol1 = 100 * sqrt(harmonic_pw) / volt1_rms;

% voltage_2
fundamental_amplitude_vol2 = Y_vol2_h(fundamental_harmonic_indx);
volt2_rms = fundamental_amplitude_vol2/sqrt(2);
harmonic_pw = sum((Y_vol2_h / sqrt(2)).^2) - (volt2_rms)^2;
THD_vol2 = 100 * sqrt(harmonic_pw) / volt2_rms;

% voltage_3
fundamental_amplitude_vol3 = Y_vol3_h(fundamental_harmonic_indx);
volt3_rms = fundamental_amplitude_vol3/sqrt(2);
harmonic_pw = sum((Y_vol3_h / sqrt(2)).^2) - (volt3_rms)^2;
THD_vol3 = 100 * sqrt(harmonic_pw) / volt3_rms;


% ..................................DPF....................................
% Calculating the angle of the fundamental harmonic of each signal
angl_curr = angle(Y_curr(fundamental_harmonic_indx));
angl_vol1 = angle(Y_vol1(fundamental_harmonic_indx));
angl_vol2 = angle(Y_vol2(fundamental_harmonic_indx));
angl_vol3 = angle(Y_vol3(fundamental_harmonic_indx));

% Calculating the DPF values
DPF_1 = cos(angl_vol1 - angl_curr);
DPF_2 = cos(angl_vol2 - angl_curr);
DPF_3 = cos(angl_vol3 - angl_curr);


%..............................Instanteneous Power.........................
% Calculating the rms values p(t)=v(t)*i(t)
IP1 = voltage1 .* current;
IP2 = voltage2 .* current;
IP3 = voltage3 .* current;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f5 = figure('Renderer', 'painters','Name','Instanteneous power','NumberTitle','off', 'Position', [10 10 1300 600]);

% Plotting the instantaneous power for voltage_1
subplot(3,1,1)
plot(time,IP1)
title('instantaneous power for voltage\_1')
xlabel('Time (sec)')
ylabel('Power (Watt)')

% Plotting the instantaneous power for voltage_2
subplot(3,1,2)
plot(time,IP2)
title('instantaneous power for voltage\_2')
xlabel('Time (sec)')
ylabel('Power (Watt)')

% Plotting the instantaneous power for voltage_3
subplot(3,1,3)
plot(time,IP3)
title('instantaneous power for voltage\_3')
xlabel('Time (sec)')
ylabel('Power (Watt)')

% Saving the figure
print -depsc Q4_voltages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%..............................Apparent Power..............................
% Calculating the rms values
V1_rms = sqrt(sum((Y_vol1_h/sqrt(2)).^2));
V2_rms = sqrt(sum((Y_vol2_h/sqrt(2)).^2));
V3_rms = sqrt(sum((Y_vol3_h/sqrt(2)).^2));
I_rms = sqrt(sum((Y_curr_h/sqrt(2)).^2));

% Calculating the S values
S1 = V1_rms * I_rms;
S2 = V2_rms * I_rms;
S3 = V3_rms * I_rms;

%...............................Active Power...............................
% Calculating the P values
P1 = sum(IP1)./L;
P2 = sum(IP2)./L;
P3 = sum(IP3)./L;

% ...................................PF....................................
% Calculating the PF values
PF1 = P1 / S1;
PF2 = P2 / S2;
PF3 = P3 / S3;


%% Question 5

% Extracting the time index of t = 0.1s (Period of signal)
t_period_indx = find(time<=0.1, 1, 'last');

% Extracting the voltage and time signals for only one Period
time_period = time(1:t_period_indx);
voltage2_period = voltage2(1:t_period_indx);

% Calculating the rms value theoreticaly 
V2_rms_theoretical = sqrt(1/length(time_period) * sum(voltage2_period.^2));



% Extracting the index of 700Hz
sample_indx = find(f==700);
 
% Calculating the rms value of Voltage2 througth the values until 700Hz
V2_rms_experimental = sqrt(sum((Y_vol2_h(1:sample_indx)/sqrt(2)).^2));;