function phase3(a_in, L_in, R, t, dt, t_period, omega, Vrms)

V_ab = Vrms*sqrt(2)*sqrt(3) * sin(omega*t);
V_bc = Vrms*sqrt(2)*sqrt(3) * sin(omega*t + deg2rad(120));
V_ca = Vrms*sqrt(2)*sqrt(3) * sin(omega*t + deg2rad(240));

plot(t_period, V_ab(1:length(t_period)),Color='blue')
hold on;
plot(t_period, V_bc(1:length(t_period)),Color='red')
hold on;
plot(t_period, V_ca(1:length(t_period)),Color='green')

hold on;
plot(t_period, -V_ab(1:length(t_period)),Color='blue',LineStyle="--")
hold on;
plot(t_period, -V_bc(1:length(t_period)),Color='red',LineStyle="--")
hold on;
plot(t_period, -V_ca(1:length(t_period)),Color='green',LineStyle="--")



period=1/50;

set(gca,'XTick',0:period/4:2*period) 
set(gca,'XTickLabel',{'0','\pi/2','\pi','3*\pi/2','2*\pi','5*\pi/2','3*\pi','7*\pi/2','4*\pi'})



    for L_i = [L_in(1) L_in(2)]
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


            % Calculating the system parameters
            A=-R/L_i;
            B=1/L_i;
            C=1;
            D=0;

            % Creating the system
            sys = ss(A, B, C, D); 

            % Creating the discrete system
            sysd = c2d(sys, dt);

            A_d = sysd.A;
            B_d = sysd.B;
            C_d = sysd.C;
            D_d = sysd.D;


        for i = 1:length(t)
            % Calculating the current phase value restricted in [0 2*pi]
            phase = mod(omega*t(i), 2*pi); 




            if(a_in <= phase && phase < a_in + pi/3)                      % 67 <= ω < 127
                V_out(i) = (-1 * (phase < pi))*2*V_ca(i) + V_ca(i);
            elseif(phase >= a_in + pi/3 && phase < a_in + 2*pi/3)           % 127 <= ω < 187
                V_out(i) = V_ab(i);
            elseif(phase >= a_in + 2*pi/3 && phase < a_in + pi)             % 187 <= ω < 247
                V_out(i) = -V_bc(i);
            elseif(phase >= a_in + pi && phase < a_in + 4*pi/3)             % 247 <= ω < 307
                V_out(i) = V_ca(i);                
            elseif(phase >= a_in + 4*pi/3 && phase < a_in + 5*pi/3) % 307 <= 367
                V_out(i) = -V_ab(i);         
            else
                V_out(i) = V_bc(i);
            end
        end  

        hold on;
        plot(t_period, V_out(1:length(t_period)),Color='black',LineStyle="-")        


        
end