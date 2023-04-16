clc;
clear;
close all;

%% Exercise data

% 
Vrms = 230;
f = 50;
R = 2.5;

Time = 10;
dt = 0.00002;

% Creating a vector for represent the total time with time step dt
t_total = 0:dt:Time;

% Creating a vector for represent the time of two periods with time step dt
t_periods = 0:dt:2/f;
% Creating a vector for represent the time until 0.4 secs with time step
% dt for representing the output current
t_I = 0:dt:0.4;

% Calculating the value of omega
omega = 2*pi*f;


%% Question 2

% Creating the input voltage for the single phase
V_i = Vrms*sqrt(2) * sin(omega*t_total);


% a = [0 deg2rad(90)];
% L = [0.04 0.08];
% 
% inp_struct = {'V_in', V_in};
% 
% for a_i = [a(1) a(2)]       
%     for L_i = [L(1) L(2)]
%         % Initializing the needed vectors to zero
%         V_out = zeros(1, length(t));
%         I_out = zeros(1, length(t));
%         x_sys = zeros(1, length(t)+1); 
% 
%         %Initializing thyristor vectors
%         T1 = zeros(1, length(t));       
%         T2 = zeros(1, length(t));       
%         T3 = zeros(1, length(t));      
%         T4 = zeros(1, length(t));  
% 
%         % Calculating the system parameters
%         A=-R/L_i;
%         B=1/L_i;
%         C=1;
%         D=0;
% 
%         % Creating the continuous and discrete system
%         sys = ss(A, B, C, D); 
%         sysd = c2d(sys, dt);
% 
%         % Extracting the variables????
%         A_d = sysd.A;
%         B_d = sysd.B;
%         C_d = sysd.C;
%         D_d = sysd.D;
% 
%         for i = 1:length(t)
%             % Calculating the current phase value restricted in [0 2*pi]
%             phase = mod(omega*t(i), 2*pi);
% 
%             % If the phase is in the terretory [a, a+pi]
%             if(phase >= a_i && phase < a_i+ pi)
%                 % Using the positive input Voltage
%                 V_out(i) = V_in(i);
%                 T1(i)=1;
%                 T2(i)=1;
%             % Else
%             else
%                 % Using the negative input voltage
%                 V_out(i) = -V_in(i);
%                 T3(i)=1;
%                 T4(i)=1;
%             end
% 
%             % Solving the system parameters
%             x_sys(i+1) = A_d*x_sys(i) + B_d*V_out(i);
%             I_out(i)   = C_d*x_sys(i) + D_d*V_out(i);
% 
%             if(x_sys(i+1) <= 0)
%                 % Setting the needed variables to zero
%                 x_sys(i+1) = 0;
%                 V_out(i) = 0;
% 
%                 % Recalculating the current
%                 I_out(i) = C_d*x_sys(i) + D_d*V_out(i);
%             end            
%         end
%         Thyristors_currs = {T1.*I_out T2.*I_out T3.*I_out T4.*I_out T4.*I_out T4.*I_out};
%         plotter(inp_struct, V_out, Thyristors_currs, t_period, a, L)
%     end
% end

%% Question 3

% Creating the three voltage inputs 
V_ab = Vrms*sqrt(2)*sqrt(3) * sin(omega*t_total);
V_bc = Vrms*sqrt(2)*sqrt(3) * sin(omega*t_total + deg2rad(120));
V_ca = Vrms*sqrt(2)*sqrt(3) * sin(omega*t_total + deg2rad(240));

% THIS IS FOR PRINTING, IGNORE
inp_struct = {'V_ab', V_ab, 'V_bc', V_bc, 'V_ac', V_ca};

% Creating tow vectors containing the values of a and L
a = [0 deg2rad(67)];
L = [0.04 0.08];

% For each a
for a_i = [a(1) a(2)]
    % for each L
    for L_i = [L(1) L(2)]
        % Initializing the needed vectors to zero
        V_out = zeros(1, length(t_total));
        I_out = zeros(1, length(t_total));
        x_sys = zeros(1, length(t_total)+1);

        %Initializing thyristor vectors
        T1 = zeros(1, length(t_total));       
        T2 = zeros(1, length(t_total));       
        T3 = zeros(1, length(t_total));      
        T4 = zeros(1, length(t_total));  
        T5 = zeros(1, length(t_total));  
        T6 = zeros(1, length(t_total));  
        
        % Calculating the system variables
        A = -R/L_i;
        B = 1/L_i;
        C = 1;
        D = 0;

        % Creating the continuous and discrete system
        sys = ss(A, B, C, D); 
        sysd = c2d(sys, dt);
        
        % Extracting the system variables
        A_d = sysd.A;   
        B_d = sysd.B;
        C_d = sysd.C;
        D_d = sysd.D;

        % For each time step
        for i = 1:length(t_total)
            % Calculating the current phase value restricted in [0 2*pi]
            phase = mod(omega*t_total(i), 2*pi);

            % According to the value of a and the value of phase the output
            % voltgage signal is created according to the thyristors that 
            % are on

            if(a_i < phase && phase <= a_i + pi/3)
                V_out(i) = -V_ca(i);
                T1(i) = 1;
                T2(i) = 1;
            elseif(a_i + pi/3 < phase && phase <= a_i + 2*pi/3)
                V_out(i) = V_ab(i);
                T1(i) = 1;
                T6(i) = 1;
            elseif(a_i + 2*pi/3 < phase && phase <= a_i + pi)
                V_out(i) = -V_bc(i);  
                T5(i) = 1;
                T6(i) = 1;
            elseif(a_i + pi <= phase && phase <= a_i + 4*pi/3)
                V_out(i) = V_ca(i);
                T4(i) = 1;
                T5(i) = 1;
            elseif(a_i + 4*pi/3 <= phase  && phase <= a_i + 5*pi/3 || a_i > pi/3 && phase <= mod(a_i + 5*pi/3, 2*pi))
                V_out(i) = -V_ab(i);
                T3(i) = 1;
                T4(i) = 1;
            else
                V_out(i) = V_bc(i);
                T2(i) = 1;
                T3(i) = 1;
            end

            % Solving the system parameters
            x_sys(i+1)  = A_d*x_sys(i)  + B_d*V_out(i);
            I_out(i)  = C_d*x_sys(i)  + D_d*V_out(i);
                
            if(x_sys(i+1) <= 0)
                % Setting the needed variables to zero
                x_sys(i+1) = 0;
                V_out(i) = 0;  
                
                % Recalculating the current
                I_out(i) = C_d*x_sys(i) + D_d*V_out(i);
            end        
        end
    % Plotting the signals
    Thyristors = {T1.*I_out T2.*I_out T3.*I_out T4.*I_out T4.*I_out T4.*I_out};
    plotter(inp_struct, V_out, I_out, Thyristors, t_I, t_periods, a_i, L_i)
    end
end

%% Question 4



