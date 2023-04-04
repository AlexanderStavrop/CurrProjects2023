clear all;
close all;
clc;

%% --------------/1-PHASE/--------------
% clear all;
% close all;
% clc;
%------------/2/------------
Vrms_1ph = 230;     % V_rms
f = 50;             % Frequency in Hz
Ts = 0.00002;       % Sampling step
T = 10;             % Simulation time
T1 =0.02           % 1-Period time
n = 500;            % Number of periods
t = 0:Ts:n*T1;      % Time vector
w = 2*pi*f;         % Frequency in rad/s
deg = w*t*180/pi;  % Angle in rad

a=[0 90];           % Vector of angle a
R=2.5;              % Resistance
L=[0.04 0.08];      % Coil vector

Vm = Vrms_1ph*sqrt(2);      % Voltage amplitude
V1ph = Vm*sin(w*t);    % 1-phase voltage

figure;
plot(deg,V1ph,'LineWidth',1);
hold on;
plot(deg,-V1ph,'--r','LineWidth',1);
axis([0 360 -400 400]);
xticks(0:30:360);
grid on;
xlabel('Angle(degrees)');
ylabel('Voltage(V)');
title('1-phase input voltage');

% Initialization of helper matrices
I_L_1phase = zeros(1,size(t,2));
V_L_1phase = zeros(1,size(t,2));
I_thyr_12 = zeros(1,size(t,2));
% I_thyr_34 = zeros(1,size(t,2));


% for k=1:2
%     for u=1:2
        
%         % Parameters of the system in space state
%         A_1phase = -R/L(u);
%         B_1phase = 1/L(u);
%         C_1phase = 1;
%         D_1phase=0;
        
%         % System in space state
%         sys_1phase = ss(A_1phase,B_1phase,C_1phase,D_1phase);
%         % Discrete system
%         sys_discrete = c2d(sys_1phase, Ts);
        
%         j=1;
        
%         % Computation of load's current and voltage and each thyristor's
%         % current
%         for i=0:Ts:n*T1
            
%             if j==1
%                 I_L_1phase(j) = 0;
%             end
            
%             omega_t = w*t(j)*180/pi;    % Angle in degrees
            
%             % Check with the angle and compute the voltage in the load
%             if mod(omega_t,2*180) < a(k)    % Angle < a
%                 V_L_1phase(j) = -V1ph(j);   % V_load = -Vin
%             elseif a(k) <= mod(omega_t,2*180) && mod(omega_t,2*180) < 180+a(k)  % a < Angle < 180+a
%                 V_L_1phase(j) = V1ph(j);    % V_load = Vin
%             else
%                 V_L_1phase(j) = -V1ph(j);   % V_load = -Vin
%             end
%             % Compute load's current
%             I_L_1phase(j+1) = sys_discrete.A*I_L_1phase(j)+sys_discrete.B*V_L_1phase(j);
            
%             % If current is negative make it 0.
%             % And the voltage as well
%             if I_L_1phase(j+1)<=0
%                 I_L_1phase(j+1) = 0;
%                 V_L_1phase(j) = 0;
%             end
            
%             % Phase 1 -> Thyr_1,2 | Phase 2 -> Thyr_3,4
%             if (mod(omega_t,2*180) >= a(k)) && (mod(omega_t,2*180) < 180 + a(k)) % Angle -> [a, 180+a)
%                 I_thyr_12(j) =  I_L_1phase(j);
%                 I_thyr_34(j) = 0;
%             else
%                 I_thyr_34(j) =  I_L_1phase(j);
%                 I_thyr_12(j) = 0;
%             end
            
%             j=j+1;
            
%         end
        
%         figure;
%         plot(t,I_L_1phase(1:end-1),'LineWidth',1);
%         grid on;
%         xlabel('Time(sec)');
%         ylabel('I_L(A)');
%         title(sprintf('1-phase load Current when a= %d, L=%.2f',a(k),L(u)));
%         axis([0.86 0.9, 0 90]);
        
%         figure;
%         subplot(1,2,1);
%         plot(t,I_thyr_12,'LineWidth',1);
%         grid on;
%         xlabel('Time(sec)');
%         ylabel('I_{T1,T2}(A)');
%         title('Thyristor(1,2)');
%         axis([0.86 0.88, 0 90]);
        
%         subplot(1,2,2);
%         plot(t,I_thyr_34,'LineWidth',1);
%         grid on;
%         xlabel('Time(sec)');
%         ylabel('I_{T3,T4}(A)');
%         title('Thyristor(3,4)');
%         axis([0.86 0.88, 0 90]);
%         suptitle(sprintf('Current in thyristors when a= %d, L=%.2f',a(k),L(u)));
        
        
%     end
    
%     figure;
%     plot(t,V_L_1phase,'LineWidth',1);
%     grid on;
%     xlabel('Time(sec)');
%     ylabel('V_L(V)');
%     title(sprintf('1-phase load Voltage when a=%d, L=0.04,0.08',a(k)));
%     axis([0.8 0.84, -350 350]);
% end






% %% --------------/3-PHASE/--------------

% %------------/3/------------

% % clear all;
% % close all;
% % clc;

% Vrms_1ph = 230;     % V_rms
% f = 50;             % Frequency in Hz
% Ts = 0.00002;       % Sampling step
% T = 10;             % Simulation time
% T1 =0.02;           % 1-Period time
% n = 500;            % Number of periods
% t = 0:Ts:n*T1;      % Time vector
% w = 2*pi*f;         % Frequency in rad/s
% deg = w*t*180/pi;   % Angle in rad

% a=[0 67];           % Vector of angle a
% R=2.5;              % Resistance
% L=[0.04 0.08];      % Coil vector

% Vm = Vrms_1ph*sqrt(2);      % Voltage amplitude
% V1ph = Vm*sin(2*pi*f*t);    % 1-phase voltage

% Van = Vm*sin(w*t - pi/6);              % Phase Voltage (Phase a)
% Vbn = Vm*sin(w*t - 2*pi/3 - pi/6);     % Phase Voltage (Phase b)
% Vcn = Vm*sin(w*t - 4*pi/3 - pi/6);     % Phase Voltage (Phase c)

% Vab = Van - Vbn;    % Voltage across phases a,b
% Vbc = Vbn - Vcn;    % Voltage across phases b,c
% Vca = Vcn - Van;    % Voltage across phases c,a
% Vba = -Vab;         % Reverse Voltage across phases a,b
% Vcb = -Vbc;         % Reverse Voltage across phases b,c
% Vac = -Vca;         % Reverse Voltage across phases c,a

% figure;
% hold on;
% plot(deg,Vab,'b','LineWidth',1);
% plot(deg,Vbc,'r','LineWidth',1);
% plot(deg,Vca,'g','LineWidth',1);
% plot(deg,Vba,'--b','LineWidth',1);
% plot(deg,Vcb,'--r','LineWidth',1);
% plot(deg,Vac,'--g','LineWidth',1);
% grid on;
% axis([0 400 -600 600]);
% xticks(0:30:360);
% xlabel('Angle(degrees)');
% ylabel('Polar voltages(V)');
% legend('Vab','Vbc','Vca','Vba','Vcb','Vac');
% title('3-phase input voltage');

% I_L_3phase = zeros(1,size(t,2));    %Initialize load's current vector
% V_L_3phase = zeros(1,size(t,2));    %Initialize load's voltage vector

% %Initialize current of thyristors vectors
% I_thyr_1 = zeros(1,size(t,2));
% I_thyr_2 = zeros(1,size(t,2));
% I_thyr_3 = zeros(1,size(t,2));
% I_thyr_4 = zeros(1,size(t,2));
% I_thyr_5 = zeros(1,size(t,2));
% I_thyr_6 = zeros(1,size(t,2));

% for k=1:2
%     for u=1:2
        
%         % Parameters of the system in space state
%         A_3phase = -R/L(u);
%         B_3phase = 1/L(u);
%         C_3phase = 1;
%         D_3phase=0;
        
%         % System in space state
%         sys_3phase = ss(A_3phase,B_3phase,C_3phase,D_3phase);
%         % Discrete system
%         sys_discrete = c2d(sys_3phase, Ts);
        
%         % Check if angle -> [0,60]
%         if ((a(k)>=0) && (a(k)<=60))
            
%             j=1;
%             % Computation of load's current and voltage and each thyristor's
%             % current
%             for i=0:Ts:n*T1
                
%                 if j==1
%                     I_L_3phase(j) = 0;
%                 end
                
%                 omega_t = w*t(j)*180/pi;    % Angle in degrees
                
%                 % Check with the angle and compute the voltage in the load
%                 % and each thyristor's current
%                 if (mod(omega_t, 2*180) < a(k)) && (mod(omega_t, 2*180) > 0) % Angle -> (0, a)
%                     V_L_3phase(j) = Vca(j);
%                     I_thyr_1(j) = 0;
%                     I_thyr_2(j) = 0;
%                     I_thyr_3(j) = 0;
%                     I_thyr_4(j) = I_L_3phase(j);
%                     I_thyr_5(j) = I_L_3phase(j);
%                     I_thyr_6(j) = 0;
%                 elseif (mod(omega_t, 2*180) >= a(k)) && (mod(omega_t, 2*180)<60+a(k)) % Angle -> [a,60+a)
%                     V_L_3phase(j) = Vcb(j);
%                     I_thyr_1(j) = 0;
%                     I_thyr_2(j) = 0;
%                     I_thyr_3(j) = 0;
%                     I_thyr_4(j) = 0;
%                     I_thyr_5(j) = I_L_3phase(j);
%                     I_thyr_6(j) = I_L_3phase(j);
%                 elseif (mod(omega_t, 2*180) >= 60 +a(k)) && (mod(omega_t, 2*180)<120 +a(k)) % Angle -> [a+60,a+120)
%                     V_L_3phase(j) = Vab(j);
%                     I_thyr_1(j) = I_L_3phase(j);
%                     I_thyr_2(j) = 0;
%                     I_thyr_3(j) = 0;
%                     I_thyr_4(j) = 0;
%                     I_thyr_5(j) = 0;
%                     I_thyr_6(j) = I_L_3phase(j);
%                 elseif (mod(omega_t, 2*180) >= 120 +a(k)) && (mod(omega_t, 2*180)<180 +a(k))     % Angle -> [a+120,a+180)
%                     V_L_3phase(j) = Vac(j);
%                     I_thyr_1(j) = I_L_3phase(j);
%                     I_thyr_2(j) = I_L_3phase(j);
%                     I_thyr_3(j) = 0;
%                     I_thyr_4(j) = 0;
%                     I_thyr_5(j) = 0;
%                     I_thyr_6(j) = 0;
%                 elseif (mod(omega_t, 2*180) >= 180 +a(k)) && (mod(omega_t, 2*180)<240 +a(k))    % Angle -> [a+180,240+a)
%                     V_L_3phase(j) = Vbc(j);
%                     I_thyr_1(j) = 0;
%                     I_thyr_2(j) = I_L_3phase(j);
%                     I_thyr_3(j) = I_L_3phase(j);
%                     I_thyr_4(j) = 0;
%                     I_thyr_5(j) = 0;
%                     I_thyr_6(j) = 0;
%                 elseif (mod(omega_t, 2*180) >= 240 +a(k)) && (mod(omega_t, 2*180)<300 +a(k))    % Angle -> [a+240,300+a)
%                     V_L_3phase(j) = Vba(j);
%                     I_thyr_1(j) = 0;
%                     I_thyr_2(j) = 0;
%                     I_thyr_3(j) = I_L_3phase(j);
%                     I_thyr_4(j) = I_L_3phase(j);
%                     I_thyr_5(j) = 0;
%                     I_thyr_6(j) = 0;
%                 elseif (mod(omega_t, 2*180) >= 300 +a(k)) && (mod(omega_t, 2*180)<=360)     % Angle -> [a+240,360)
%                     V_L_3phase(j) = Vca(j);
%                     I_thyr_1(j) = 0;
%                     I_thyr_2(j) = 0;
%                     I_thyr_3(j) = 0;
%                     I_thyr_4(j) = I_L_3phase(j);
%                     I_thyr_5(j) = I_L_3phase(j);
%                     I_thyr_6(j) = 0;
%                 end
                
%                 % Compute load's current
%                 I_L_3phase(j+1) = sys_discrete.A*I_L_3phase(j)+sys_discrete.B*V_L_3phase(j);
                
%                 % If current is negative make it 0.
%                 % And the voltage as well
%                 if I_L_3phase(j+1)<=0
%                     I_L_3phase(j+1) = 0;
%                     V_L_3phase(j) = 0;
%                 end
                
                
%                 j=j+1;
%             end
            
%             % Check if angle -> (60,90]
%         elseif (a(k)>60 && a(k)<=90)
            
%             j=1;
%             % Computation of load's current and voltage and each thyristor's
%             % current
%             for i=0:Ts:n*T1
                
%                 if j==1
%                     I_L_3phase(j) = 0;
%                 end
                
%                 omega_t = w*t(j)*180/pi;    %Angle in degrees
                
%                 % Check with the angle and compute the voltage in the load
%                 if (mod(omega_t, 2*180) < a(k)-60) && (mod(omega_t, 2*180) >= 0)    % Angle -> [0, a-60)
%                     V_L_3phase(j) = Vba(j);
%                     I_thyr_1(j) = 0;
%                     I_thyr_2(j) = 0;
%                     I_thyr_3(j) = I_L_3phase(j);
%                     I_thyr_4(j) = I_L_3phase(j);
%                     I_thyr_5(j) = 0;
%                     I_thyr_6(j) = 0;
%                 elseif (mod(omega_t, 2*180) >= a(k)-60) && (mod(omega_t, 2*180)<a(k))   % Angle -> [a-60,a)
%                     V_L_3phase(j) = Vca(j);
%                     I_thyr_1(j) = 0;
%                     I_thyr_2(j) = 0;
%                     I_thyr_3(j) = 0;
%                     I_thyr_4(j) = I_L_3phase(j);
%                     I_thyr_5(j) = I_L_3phase(j);
%                     I_thyr_6(j) = 0;
%                 elseif (mod(omega_t, 2*180) >= a(k)) && (mod(omega_t, 2*180)<60+a(k))   % Angle -> [a,a+60)
%                     V_L_3phase(j) = Vcb(j);
%                     I_thyr_1(j) = 0;
%                     I_thyr_2(j) = 0;
%                     I_thyr_3(j) = 0;
%                     I_thyr_4(j) = 0;
%                     I_thyr_5(j) = I_L_3phase(j);
%                     I_thyr_6(j) = I_L_3phase(j);
%                 elseif (mod(omega_t, 2*180) >= 60 +a(k)) && (mod(omega_t, 2*180)<120 +a(k)) % Angle -> [a+60,a+120)
%                     V_L_3phase(j) = Vab(j);
%                     I_thyr_1(j) = I_L_3phase(j);
%                     I_thyr_2(j) = 0;
%                     I_thyr_3(j) = 0;
%                     I_thyr_4(j) = 0;
%                     I_thyr_5(j) = 0;
%                     I_thyr_6(j) = I_L_3phase(j);
%                 elseif (mod(omega_t, 2*180) >= 120 +a(k)) && (mod(omega_t, 2*180)<180 +a(k))    % Angle -> [a+120,a+180)
%                     V_L_3phase(j) = Vac(j);
%                     I_thyr_1(j) = I_L_3phase(j);
%                     I_thyr_2(j) = I_L_3phase(j);
%                     I_thyr_3(j) = 0;
%                     I_thyr_4(j) = 0;
%                     I_thyr_5(j) = 0;
%                     I_thyr_6(j) = 0;
%                 elseif (mod(omega_t, 2*180) >= 180 +a(k)) && (mod(omega_t, 2*180)<240 +a(k))    % Angle -> [a+180,a+240)
%                     V_L_3phase(j) = Vbc(j);
%                     I_thyr_1(j) = 0;
%                     I_thyr_2(j) = I_L_3phase(j);
%                     I_thyr_3(j) = I_L_3phase(j);
%                     I_thyr_4(j) = 0;
%                     I_thyr_5(j) = 0;
%                     I_thyr_6(j) = 0;
%                 elseif (mod(omega_t, 2*180) >= 240 +a(k)) && (mod(omega_t, 2*180)<300+a(k))     % Angle -> [a+240,300+a)
%                     V_L_3phase(j) = Vba(j);
%                     I_thyr_1(j) = 0;
%                     I_thyr_2(j) = 0;
%                     I_thyr_3(j) = I_L_3phase(j);
%                     I_thyr_4(j) = I_L_3phase(j);
%                     I_thyr_5(j) = 0;
%                     I_thyr_6(j) = 0;
%                 end
                
%                 % Compute load's current
%                 I_L_3phase(j+1) = sys_discrete.A*I_L_3phase(j)+sys_discrete.B*V_L_3phase(j);
                
%                 % If current is negative make it 0.
%                 % And the voltage as well
%                 if I_L_3phase(j+1)<0
%                     I_L_3phase(j+1) = 0;
%                     V_L_3phase(j) = 0;
%                 end
                
%                 j=j+1;
                
%             end
            
%             % Check if angle -> (90,360]
%         else
%             fprintf('Angle a must be between 0 and 90');
            
%         end
        
%         figure;
%         plot(t,I_L_3phase(1:end-1),'LineWidth',1);
%         grid on;
%         xlabel('Time(sec)');
%         ylabel('I_L(A)');
%         title(sprintf('3-phase load Current when a= %d, L=%.2f',a(k),L(u)));
%         axis([0.86 0.88 ,0 250]);  
        
        
%         figure;
%         subplot(3,2,1);
%         plot(t,I_thyr_1,'LineWidth',1);
%         grid on;
%         xlabel('Time(sec)');
%         ylabel('I_{T1}(A)');
%         title('Thyristor(1)');        
%         axis([0.86 0.88 ,0 250]);  
        
%         subplot(3,2,2)
%         plot(t,I_thyr_2,'LineWidth',1);
%         grid on;
%         xlabel('Time(sec)');
%         ylabel('I_{T2}(A)');
%         title('Thyristor(2)');
%         axis([0.86 0.88 ,0 250]);  
        
%         subplot(3,2,3)
%         plot(t,I_thyr_3,'LineWidth',1);
%         grid on;
%         xlabel('Time(sec)');
%         ylabel('I_{T3}(A)');
%         title('Thyristor(3)');
%         axis([0.86 0.88 ,0 250]);  
        
%         subplot(3,2,4)
%         plot(t,I_thyr_4,'LineWidth',1);
%         grid on;
%         xlabel('Time(sec)');
%         ylabel('I_{T4}(A)');
%         title('Thyristor(4)');
%         axis([0.86 0.88 ,0 250]);  
        
%         subplot(3,2,5)
%         plot(t,I_thyr_5,'LineWidth',1);
%         grid on;
%         xlabel('Time(sec)');
%         ylabel('I_{T5}(A)');
%         title('Thyristor(5)');
%         axis([0.86 0.88 ,0 250]);  
        
%         subplot(3,2,6)
%         plot(t,I_thyr_6,'LineWidth',1);
%         grid on;
%         xlabel('Time(sec)');
%         ylabel('I_{T6}(A)');
%         title('Thyristor(6)');
%         axis([0.86 0.88 ,0 250]);  
%         suptitle (sprintf('Current in thyristors when a= %d, L=%.2f',a(k),L(u)));
              
        
%     end
    
%     figure;
%     plot(t,V_L_3phase,'LineWidth',1);
%     grid on;
%     xlabel('Time(sec)');
%     ylabel('V_L(V)');
%     title(sprintf('3-phase load Voltage when a=%d, L=0.04,0.08',a(k)));
%     axis([0.8 0.84, -100 570]);
    
% end


% %% ------------/4.Controller/------------
% % clear all;
% % close all;
% % clc;

% Vrms_1ph = 230;         % Initial V_rms value
% f = 50;                 % Frequency in Hz
% Ts = 0.00002;           % Sampling step
% T = 10;                 % Simulation time
% T1 =0.02;               % 1-Period time
% n = 500;                % Number of periods
% t = 0:Ts:n*T1;          % Time vector
% t1 = 0:Ts:n*T1/2-Ts;    % Time vector for the first 5 secs
% w = 2*pi*f;             % Frequency in rad/s
% deg = w*t*180/pi;       % Angle in rad
% Vm = Vrms_1ph*sqrt(2);  % Initial voltage amplitude


% Van = Vm*sin(w*t - pi/6);              % Phase Voltage (Phase a)
% Vbn = Vm*sin(w*t - 2*pi/3 - pi/6);     % Phase Voltage (Phase b)
% Vcn = Vm*sin(w*t - 4*pi/3 - pi/6);     % Phase Voltage (Phase c)

% Vab = Van - Vbn;    % Voltage across phases a,b
% Vbc = Vbn - Vcn;    % Voltage across phases b,c
% Vca = Vcn - Van;    % Voltage across phases c,a
% Vba = -Vab;         % Reverse Voltage across phases a,b
% Vcb = -Vbc;         % Reverse Voltage across phases b,c
% Vac = -Vca;         % Reverse Voltage across phases c,a

% I_ref = 75; % Reference value of current
% R= 2.5;     % Resistance
% L= 0.04;    % Coil vector

% % Parameters of the system in space state
% A = -R/L;
% B = 1/L;
% C = 1;
% D = 0;

% % System in space state
% sys_3phase = ss(A,B,C,D);
% % Discrete system
% sys_discrete = c2d(sys_3phase, Ts);

% % Parameters of the controller's system in space state
% A_ctrl = 0;
% B_ctrl = 1;
% C_ctrl = 6;
% D_ctrl = 1;

% % Controller's system in space state
% sys_controller = ss(A_ctrl,B_ctrl,C_ctrl,D_ctrl);
% % Discrete system
% sys_discrete_ctrl = c2d(sys_controller, Ts);

% %Initialize helping vectors
% err = zeros(1,size(t,2));
% I_L_3phase = zeros(1,size(t,2)+1);
% V_L_3phase = zeros(1,size(t,2));
% x = zeros(1,size(t,2));
% a_4 = zeros(1,size(t,2));
% I_L_mean = zeros(1,size(t,2));

% j=1;

% for i = t1-Ts
    
%     omega_t = w*t1(j)*180/pi;   %Angle in degrees
    
%     if j==1
%         if (mod(omega_t, 2*180) < a_4(j)) && (mod(omega_t, 2*180) > 0)  % Angle ->(0,a)
%             V_L_3phase(j) = Vca(j);
%         elseif (mod(omega_t, 2*180) >= a_4(j)) && (mod(omega_t, 2*180)<60+ a_4(j))  % Angle -> [a,60+a)
%             V_L_3phase(j) = Vcb(j);
%         elseif (mod(omega_t, 2*180) >= 60 + a_4(j)) && (mod(omega_t, 2*180)<120 + a_4(j))   % Angle -> [60+a,120+a)
%             V_L_3phase(j) = Vab(j);
%         elseif (mod(omega_t, 2*180) >= 120 + a_4(j)) && (mod(omega_t, 2*180)<180 + a_4(j))  % Angle -> [120+a,180+a)
%             V_L_3phase(j) = Vac(j);
%         elseif (mod(omega_t, 2*180) >= 180 + a_4(j)) && (mod(omega_t, 2*180)<240 + a_4(j))  % Angle -> [180+a,240+a)
%             V_L_3phase(j) = Vbc(j);
%         elseif (mod(omega_t, 2*180) >= 240 +a_4(j)) && (mod(omega_t, 2*180)<300 + a_4(j))   % Angle -> [240+a,300+a)
%             V_L_3phase(j) = Vba(j);
%         elseif (mod(omega_t, 2*180) >= 300 + a_4(j)) && (mod(omega_t, 2*180)<=360)  % Angle -> [300+a,360]
%             V_L_3phase(j) = Vca(j);
%         end
        
%     else
        
%         if (mod(omega_t, 2*180) < a_4(j-1)) && (mod(omega_t, 2*180) > 0)    % Angle -> (0,a)
%             V_L_3phase(j) = Vca(j);
%         elseif (mod(omega_t, 2*180) >= a_4(j-1)) && (mod(omega_t, 2*180)<60+ a_4(j-1))  % Angle -> [a,60+a)
%             V_L_3phase(j) = Vcb(j);
%         elseif (mod(omega_t, 2*180) >= 60 + a_4(j-1)) && (mod(omega_t, 2*180)<120 + a_4(j-1))   % Angle -> [60+a, 120+a)
%             V_L_3phase(j) = Vab(j);
%         elseif (mod(omega_t, 2*180) >= 120 + a_4(j-1)) && (mod(omega_t, 2*180)<180 + a_4(j-1))  % Angle -> [120+a, 180+a)
%             V_L_3phase(j) = Vac(j);
%         elseif (mod(omega_t, 2*180) >= 180 + a_4(j-1)) && (mod(omega_t, 2*180)<240 + a_4(j-1))  % Angle -> [180+a, 240+a)
%             V_L_3phase(j) = Vbc(j);
%         elseif (mod(omega_t, 2*180) >= 240 +a_4(j-1)) && (mod(omega_t, 2*180)<300 + a_4(j-1))   % Angle -> [240+a, 300+a)
%             V_L_3phase(j) = Vba(j);
%         elseif (mod(omega_t, 2*180) >= 300 + a_4(j-1)) && (mod(omega_t, 2*180)<=360)    % Angle -> [300+a, 360]
%             V_L_3phase(j) = Vca(j);
%         end
        
%     end
%     % Compute load's current
%     I_L_3phase(j+1) = sys_discrete.A*I_L_3phase(j)+sys_discrete.B*V_L_3phase(j);
    
%     % Find the appropriate counters for the computation of I_L_mean
%     if(j<=T1/Ts)
%         min_pos=1;
%     else
%         min_pos=j-(T1/Ts-1);
%     end
    
%     I_L_mean(j) = mean(I_L_3phase(min_pos:j));  % Mean value of load's current
%     err(j) = I_L_mean(j) - I_ref;   % Calculation of error
    
%     x(j+1) = sys_discrete_ctrl.A*x(j) + sys_discrete_ctrl.B*err(j);
%     a_4(j) = sys_discrete_ctrl.C*x(j) + sys_discrete_ctrl.D*err(j);
    
%     % Check if angle belongs to [0,90]
%     if a_4(j) > 90
%         a_4(j) = 90;
%     elseif a_4(j) < 0
%         a_4(j) = 0;
%     end
    
    
%     j=j+1;
    
% end

% Vrms_1ph = 190;         % Final V_rms value
% t2 = n*T1/2:Ts:n*T1;    % Time vector for the last 5 secs
% Vm = Vrms_1ph*sqrt(2);  % Final voltage amplitude

% Van = Vm*sin(w*t - pi/6);              % Phase Voltage (Phase a)
% Vbn = Vm*sin(w*t - 2*pi/3 - pi/6);     % Phase Voltage (Phase b)
% Vcn = Vm*sin(w*t - 4*pi/3 - pi/6);     % Phase Voltage (Phase c)

% Vab = Van - Vbn;    % Voltage across phases a,b
% Vbc = Vbn - Vcn;    % Voltage across phases b,c
% Vca = Vcn - Van;    % Voltage across phases c,a
% Vba = -Vab;         % Reverse Voltage across phases a,b
% Vcb = -Vbc;         % Reverse Voltage across phases b,c
% Vac = -Vca;         % Reverse Voltage across phases c,a

% for i=t2-Ts
    
    
%     omega_t = w*t(j)*180/pi;
%     if (mod(omega_t, 2*180) < a_4(j-1)) && (mod(omega_t, 2*180) > 0)    % Angle -> (0,a)
%         V_L_3phase(j) = Vca(j);
%     elseif (mod(omega_t, 2*180) >= a_4(j-1)) && (mod(omega_t, 2*180)<60+ a_4(j-1))  % Angle -> [a,60+a)
%         V_L_3phase(j) = Vcb(j);
%     elseif (mod(omega_t, 2*180) >= 60 + a_4(j-1)) && (mod(omega_t, 2*180)<120 + a_4(j-1))   % Angle -> [60+a, 120+a)
%         V_L_3phase(j) = Vab(j);
%     elseif (mod(omega_t, 2*180) >= 120 + a_4(j-1)) && (mod(omega_t, 2*180)<180 + a_4(j-1))  % Angle -> [120+a, 180+a)
%         V_L_3phase(j) = Vac(j);
%     elseif (mod(omega_t, 2*180) >= 180 + a_4(j-1)) && (mod(omega_t, 2*180)<240 + a_4(j-1))  % Angle -> [180+a, 240+a)
%         V_L_3phase(j) = Vbc(j);
%     elseif (mod(omega_t, 2*180) >= 240 +a_4(j-1)) && (mod(omega_t, 2*180)<300 + a_4(j-1))   % Angle -> [240+a, 300+a)
%         V_L_3phase(j) = Vba(j);
%     elseif (mod(omega_t, 2*180) >= 300 + a_4(j-1)) && (mod(omega_t, 2*180)<=360)    % Angle -> [300+a, 360]
%         V_L_3phase(j) = Vca(j);
%     end
    
%     % Compute load's current
%     I_L_3phase(j+1) = sys_discrete.A*I_L_3phase(j)+sys_discrete.B*V_L_3phase(j);
    
%     % Find the appropriate counters for the computation of I_L_mean
%     if(j<=T1/Ts)
%         min_pos=1;
%     else
%         min_pos=j-(T1/Ts-1);
%     end
    
%     I_L_mean(j) = mean(I_L_3phase(min_pos:j));  % Mean value of load's current
%     err(j) = I_L_mean(j) - I_ref;   % Calculation of error
    
%     %Calculation of the angle (PI controllers output)
%     x(j+1) = sys_discrete_ctrl.A*x(j) + sys_discrete_ctrl.B*err(j);
%     a_4(j) = sys_discrete_ctrl.C*x(j) + sys_discrete_ctrl.D*err(j);
    
%     % Check if angle belongs to [0,90]
%     if a_4(j) > 90
%         a_4(j) = 90;
%     elseif a_4(j) < 0
%         a_4(j) = 0;
%     end
    
%     j=j+1;
    
% end


% figure;
% plot(t,I_L_3phase(1:end-1),'LineWidth',1);
% grid on;
% hold on
% plot(t,75*ones(1,size(t,2)),'LineWidth',1);
% xlabel('Time(secs)');
% ylabel('I_L(A)');
% title(sprintf('3-phase load Current when R=%.2f, L=%.2f',R,L));
% legend('I_{Load}','I_{ref}');

% figure;
% plot(t,a_4,'LineWidth',1);
% grid on;
% hold on
% xlabel('Time(sec)');
% ylabel('a(deg)');
% title('Angle of thyristors');

% figure;
% plot(t,I_L_mean,'LineWidth',1);
% grid on;
% hold on;
% plot(t,75*ones(1,size(t,2)),'LineWidth',1);
% xlabel('Time(sec)');
% ylabel('I_{L,mean}(A)');
% title(sprintf('Mean Current in the load when R= %.2f, L=%.2f',R,L));
% legend('I_{Load,mean}','I_{ref}');

% figure;
% plot(t,V_L_3phase,'LineWidth',1);
% grid on;
% xlabel('Time(sec)');
% ylabel('V_L(V)');
% title(sprintf('3-phase load Voltage when R=%.2f, L=%.2f',R,L));





