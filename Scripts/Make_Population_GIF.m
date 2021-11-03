% Code for making GIF of carrier population change during a JV sweep with
% and without ion motion

%% Read in data

par_10 = pc('Input_files/ptaa_mapi_pcbm_with_vsr.csv');

%% Find eqm solution

eqm = equilibrate(par_10);

%% Measure JV curves with and wihtout ion motion

CV_el = doCV(eqm.el, 1, 0, 1.3, 0, 10e-3, 1, 241);
CV_ion = doCV(eqm.ion, 1, 0, 1.3, 0, 10e-3, 1, 241);

%% Do the plotting 
figure('WindowState', 'fullscreen')

filename = 'C:\Users\User\Documents\Work\PhD\Figures\Populations.gif';
for n = 1:241
    % Define figure and axes
    plot(CV_el.x, log10(CV_el.u(n,:,2)), 'b--')

    hold on

    % Do background colours
    % Arguments of patch give the coordinates of the corners of the polygon to
    % be shaded
    patch([0 CV_el.x(160) CV_el.x(160) 0], [18 18 0 0], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
    patch([CV_el.x(160) CV_el.x(260) CV_el.x(260) CV_el.x(160)], [18 18 0 0], 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'none')
    patch([CV_el.x(660) CV_el.x(760) CV_el.x(760) CV_el.x(660)], [18 18 0 0], 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'none')
    patch([CV_el.x(760) CV_el.x(end) CV_el.x(end) CV_el.x(760)], [18 18 0 0],  'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none')

    %Plot actual data
    plot(CV_el.x, log10(CV_el.u(n,:,2)), 'b--')
    plot(CV_el.x, log10(CV_el.u(n,:,3)), 'r--')
    plot(CV_ion.x, log10(CV_ion.u(n,:,2)),'b')
    plot(CV_ion.x, log10(CV_ion.u(n,:,3)), 'r')

    hold off

    % Labels and tick marks
    xlim([0, sum(CV_el.par.d)])
    xticks([0 100e-7 200e-7 300e-7 400e-7 500e-7])
    xticklabels({'0', '100', '200', '300', '400', '500'})
    ylim([8,18])
    yticks([8 10 12 14 16 18])
    yticklabels({'10^{8}', '10^{10}', '10^{12}', '10^{14}', '10^{16}', '10^{18}'})
    xlabel('distance (nm)')
    ylabel('carrier concentration (cm^{-3})')
    legend({'','','','','',' n, No Ions', ' p, No Ions',' n, Ions', ' p, Ions' }, 'Position', [0.75 0.8 0.1 0.1])

    %Make gif of populatin change across JV curve
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if n == 1
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end