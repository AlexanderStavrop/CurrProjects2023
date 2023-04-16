function plotter(volts, vout, iout, thyristors, t_I, t_periods, a, L)
    path = '~/Downloads/Exercise_2/Images/';
    colors = ["cyan" "red" "green" "black"];

    %%
    % ΝΑ ΦΤΙΑΞΩ ΤΟΥΣ ΑΞΟΝΕΣ ΚΑΙ ΤΑ ΟΝΟΜΑΤΑ ΣΤΑ ΡΕΥΜΑΤΑ !!!!!!!!!!!!!!!!!!!!
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig = figure('Renderer', 'painters','Name','V_in vs V_out','NumberTitle','off', 'Position', [10 10 900 540] );
    title(sprintf('V_{out} vs V_{in} (a=%s deg, L=%.2f H)', num2str(a*180/pi), L))
    set(0,'DefaultLineLineWidth',1.2)

    for i = 1:2:length(volts)
        name = volts{i};
        name = sprintf('%s{%s}', name(1:2), name(3:end));
        signal = volts{i + 1};

        plot(t_periods, signal(1:length(t_periods)), 'color', colors(ceil(i/2)), 'DisplayName',name)
        hold on;
        name = sprintf("-%s", name);
        plot(t_periods, -signal(1:length(t_periods)), 'color', colors(ceil(i/2)), 'LineStyle', '--', 'DisplayName',name)
        hold on;
    end

    plot(t_periods, vout(1:length(t_periods)), 'color', colors(4), "LineStyle"','--', 'DisplayName','V_{out}')
    legend()

    fname = sprintf('%s3_Vout_%s_%s',path, num2str(a*180/pi), sprintf('%02d', L*100));
    print(fname, '-depsc')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig = figure('Renderer', 'painters','Name','I_out','NumberTitle','off', 'Position', [10 10 900 540]);
    title(sprintf('I_{out} (a=%s deg, L=%.2f H)', num2str(a*180/pi), L))

    plot(t_I, iout(1:length(t_I)), 'DisplayName','I_{out}')
    legend()

    xstart=.4;
    xend=0.87;
    ystart=.2;
    yend=.5;
    axes('position',[xstart ystart xend-xstart yend-ystart ])
    box on

    start_indx = find(t_I<=0.28, 1, 'last');
    end_indx = find(t_I<=0.32, 1, 'last');

    plot(t_I(start_indx:end_indx), iout((start_indx:end_indx)), 'DisplayName','I_{out}')
    fname = sprintf('%s3_Iout_%s_%s',path, num2str(a*180/pi), sprintf('%02d', L*100));
    print(fname, '-depsc')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig = figure('Renderer', 'painters','Name','Thyristor_currents','NumberTitle','off', 'Position', [10 10 1300 600] );
    title(sprintf('Thyristor_currents (a=%s deg, L=%.2f H)', num2str(a*180/pi), L))
    set(0,'DefaultLineLineWidth',1.2)

    numOfThyristors = length(thyristors);
    for i = 1:length(thyristors)
        subplot(2, numOfThyristors/2, i)

        curr = thyristors{i};
        plot(t_periods,curr(1:length(t_periods)))
        xlabel('Time (sec)');
        ylabel('Cuurent Amplitude (A)');
    end

    fname = sprintf('%s3_ThI_%s_%s',path, num2str(a*180/pi), sprintf('%02d', L*100));
    print(fname, '-depsc')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%