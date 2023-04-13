clc;
clear;
close all;

%% Question 1
Vrms=230;
f=50;
dt=0.00002;
T=10;

t=0:dt:T;
t_period = 0:dt:1/f;

omega = 2*pi*f;

V_in = Vrms*sqrt(2) * sin(omega*t);
R=2.5;


%% Question 2
a = [0 deg2rad(90)];
L = [0.04 0.08];

inp_struct = {'V_in', V_in};

for a_i = [a(1) a(2)]       
    for L_i = [L(1) L(2)]
        % Initializing the needed vectors to zero
        V_out = zeros(1, length(t));
        I_out = zeros(1, length(t));
        x_sys = zeros(1, length(t)+1); 
        
        %Initializing thyristor vectors
        T1 = zeros(1, length(t));       
        T2 = zeros(1, length(t));       
        T3 = zeros(1, length(t));      
        T4 = zeros(1, length(t));  

        % Calculating the system parameters
        A=-R/L_i;
        B=1/L_i;
        C=1;
        D=0;
        
        % Creating the continuous and discrete system
        sys = ss(A, B, C, D); 
        sysd = c2d(sys, dt);
        
        % Extracting the variables????
        A_d = sysd.A;
        B_d = sysd.B;
        C_d = sysd.C;
        D_d = sysd.D;
        
        for i = 1:length(t)
            % Calculating the current phase value restricted in [0 2*pi]
            phase = mod(omega*t(i), 2*pi);
            
            % If the phase is in the terretory [a, a+pi]
            if(phase >= a_i && phase < a_i+ pi)
                % Using the positive input Voltage
                V_out(i) = V_in(i);
                T1(i)=1;
                T2(i)=1;
                T3(i)=0;
                T4(i)=0;
            % Else
            else
                % Using the negative input voltage
                V_out(i) = -V_in(i);
                T1(i)=0;
                T2(i)=0;
                T3(i)=1;
                T4(i)=1;
            end
            
            % Solving the system parameters
            x_sys(i+1) = A_d*x_sys(i) + B_d*V_out(i);
            I_out(i)   = C_d*x_sys(i) + D_d*V_out(i);

            if(x_sys(i+1) <= 0)
                % Setting the needed variables to zero
                x_sys(i+1) = 0;
                V_out(i) = 0;
                
                % Recalculating the current
                I_out(i) = C_d*x_sys(i) + D_d*V_out(i);
            end            
        end
        Thyristors_currs = {T1.*I_out T2.*I_out T3.*I_out T4.*I_out T4.*I_out T4.*I_out};
        plotter(inp_struct, V_out, Thyristors_currs, t_period, a, L)
    end
end

%% Question 3

V_ab = Vrms*sqrt(2)*sqrt(3) * sin(omega*t);
V_bc = Vrms*sqrt(2)*sqrt(3) * sin(omega*t + deg2rad(120));
V_ca = Vrms*sqrt(2)*sqrt(3) * sin(omega*t + deg2rad(240));

inp_struct = {'V_ab', V_ab, 'V_bc', V_bc, 'V_ca', V_ca};

% period=1/50;
% 
% set(gca,'XTick',0:period/4:2*period) 
% set(gca,'XTickLabel',{'0','\pi/2','\pi','3*\pi/2','2*\pi','5*\pi/2','3*\pi','7*\pi/2','4*\pi'})



a = [0 deg2rad(67)];
L = [0.04 0.08];

for a_i = [a(1) a(2)]
    for L_i = [L(1) L(2)]
        % Initializing the needed vectors to zero
        V_out = zeros(1, length(t));
        I_out = zeros(1, length(t));
        x_sys = zeros(1, length(t)+1); 

        %Initializing thyristor vectors
        T1 = zeros(1, length(t));       
        T2 = zeros(1, length(t));       
        T3 = zeros(1, length(t));      
        T4 = zeros(1, length(t));  
        T5 = zeros(1, length(t));  
        T6 = zeros(1, length(t));  
        
        A=-R/L_i;
        B=1/L_i;
        C=1;
        D=0;

        % Creating the continuous and discrete system
        sys = ss(A, B, C, D); 
        sysd = c2d(sys, dt);
        
        % Extracting the variables????
        A_d = sysd.A;
        B_d = sysd.B;
        C_d = sysd.C;
        D_d = sysd.D;

        for i = 1:length(t)
            % Calculating the current phase value restricted in [0 2*pi]
            phase = mod(omega*t(i), 2*pi); 
%             bool = (phase < pi);
%             not_bool = (phase >= pi);
%             if(a_i + not_bool*pi <= phase && phase < a_i + pi/3 + not_bool*pi)
%                 V_out(i) = -bool*2*V_ca(i) + V_ca(i);
%             
%             elseif(a_i + pi/3 + not_bool*pi <= phase && phase < a_i + 2*pi/3 + not_bool*pi)           
%                 V_out(i) = -not_bool*2*V_ab(i) + V_ab(i);
% 
%             elseif(a_i + 2*pi/3 + not_bool*pi <= phase && phase < a_i + pi + not_bool*pi)             
%                 V_out(i) = -bool*2*V_bc(i) + V_bc(i);
%             end

            if(a_i <= phase && phase < a_i + pi/3)
                V_out(i) = (-1 * (phase < pi))*2*V_ca(i) + V_ca(i);
            elseif(phase >= a_i + pi/3 && phase < a_i + 2*pi/3)
                V_out(i) = V_ab(i);
            elseif(phase >= a_i + 2*pi/3 && phase < a_i + pi)
                V_out(i) = -V_bc(i);
            elseif(phase >= a_i + pi && phase < a_i + 4*pi/3)
                V_out(i) = V_ca(i);
            elseif(phase >= a_i + 4*pi/3 && phase < a_i + 5*pi/3)
                V_out(i) = -V_ab(i);
            else
                V_out(i) = V_bc(i);
            end
        end
        
        plotter(inp_struct, V_out, '0', t_period, a, L)
    end
end