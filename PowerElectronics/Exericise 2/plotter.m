function plotter(volts, vout, thyristors, t_period, a, L)
    
    colors = ["cyan" "red" "green" "black"];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig = figure('Renderer', 'painters','Name','V_in vs V_out','NumberTitle','off', 'Position', [10 10 1300 600] );
    title(sprintf('V_{out} vs V_{in} (a=%s deg, L=%.2f H)', num2str(a*180/pi), L))
    set(0,'DefaultLineLineWidth',1.2)
    
    for i = 1:2:length(volts)
        name = volts{i};
        signal = volts{i + 1};

        plot(t_period, signal(1:length(t_period)), 'color', colors(ceil(i/2)), 'DisplayName',name)
        hold on;
        plot(t_period, -signal(1:length(t_period)), 'color', colors(ceil(i/2)), 'LineStyle', '--')
        hold on;
    end
    
    plot(t_period, vout(1:length(t_period)), 'color', colors(4), "LineStyle"','--')
    legend()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     fig = figure('Renderer', 'painters','Name','V_in vs V_out','NumberTitle','off', 'Position', [10 10 1300 600] );
%     title(sprintf('V_{out} vs V_{in} (a=%s deg, L=%.2f H)', num2str(a*180/pi), L))
%     figure()
%     set(0,'DefaultLineLineWidth',1.2)
% 
%     numOfThyristors = length(thyristors);
%     for i = 1:length(thyristors)
%         subplot(2, numOfThyristors/2, i)
%         
%         curr = thyristors{i};
%         plot(t_period,curr(1:length(t_period)))
%         xlabel('Time (sec)');
%         ylabel('Cuurent Amplitude (A)');
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%