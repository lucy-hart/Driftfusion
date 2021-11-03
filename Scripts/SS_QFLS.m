%Script for investigating amount of QFLS at SC and quenching of PL between
%OC and SC as a function of interlayer material and width
%Parameters and CV scan rate taken from 'Voltage-Dependent
%Photoluminescence and How It Correlates with the Fill Factor 
%and Open-Circuit Voltage in Perovskite Solar Cells'.

%% Read in data and create devices with HTLs of different thicknesses
par_10 = pc('Input_files/ptaa_mapi_pcbm_with_vsr.csv');

%par_10 = pc('Input_files/ptaa_mapi_pcbm.csv');
%par_10.vsr_mode=0;
%par_10 = refresh_device(par_10);

par_50 = par_10;
par_50.d(1) = 5*10^-6;

par_50 = refresh_device(par_50);

par_50_doped = par_50;
par_50_doped.E0(1) = -5.3;

par_50_doped = refresh_device(par_50_doped);

par_90 = par_10;
par_90.d(1) = 9*10^-6;
par_90 = refresh_device(par_90);

devices = {par_10, par_50, par_90, par_50_doped};

%% Find SS solutions at 1 sun intensity
eqm_solutions_dark = cell(1,4);
eqm_solutions_light = cell(1,4);
for i = 1:4
    eqm = equilibrate(devices{i});
    eqm_solutions_dark{i} = eqm;
    eqm_solutions_light{i} = changeLight(eqm.ion,1,0);
end

%% Perform CV scans
CV_solutions = cell(1,4);
for j = 1:4
    sol = eqm_solutions_dark{j}.el;
    CV_solutions{j} = doCV(sol, 1, 0, 1.3, 0, 10e-3, 1, 241);
    dfplot.JtotVapp(CV_solutions{j},0)
    hold on
end
hold off
legend({'10 nm', '50 nm', '90 nm', '50 nm, doped'}, 'Location', 'northwest')
xlim([0, 1.3])
ylim([-0.03, 0.01])

%% Break down contributions to the current
%Columns in J_values are J_gen, J_rad, J_srh, J_vsr and J_ext
J_values = zeros(241,5,4);
e = 1.60e-19;
for k = 1:4
    loss_currents = dfana.calcr(CV_solutions{k},'whole');
    J_values(:,1,k) = -e*trapz(CV_solutions{k}.par.x_sub, dfana.calcg(CV_solutions{k}).').';
    J_values(:,2,k) = e*trapz(CV_solutions{k}.par.x_sub, loss_currents.btb(:,1:end-1).').';
    J_values(:,3,k) = e*trapz(CV_solutions{k}.par.x_sub, loss_currents.srh(:,1:end-1).').';
    J_values(:,4,k) = e*trapz(CV_solutions{k}.par.x_sub, loss_currents.vsr(:,1:end-1).').';
    J_values(:,5,k) = dfana.calcJ(CV_solutions{k}).tot(:,1);
    
end

%% Plot contributons to the current
figure(2)
line_colour = {[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0 0.4470 0.7410], [0.3010 0.7450 0.9330], [0.4660 0.6740 0.1880]};
V_10 = dfana.calcVapp(CV_solutions{1});
V_90 = dfana.calcVapp(CV_solutions{3});
hold on
for n = 1:5
    plot(V_10(1:121), J_values(1:121,n,3), 'color', line_colour{n}) 
    plot(V_90(1:121), J_values(1:121,n,1), 'color', line_colour{n}, 'linestyle', '--')
end
plot(V_10(1:121), zeros(1,121), 'black', 'LineWidth', 1)
hold off
xlim([0, 1.3])
xlabel('Voltage (V)')
ylim([-0.025, 0.025])
ylabel('Current Density (Acm^{-2})')
legend({'J_{gen}', '', 'J_{rad}', '', 'J_{SRH}', '', 'J_{VSR}','', 'J_{ext}', '', ''}, 'Location', 'bestoutside')

%% Plot 'PL' results
figure(11)
for i = 1:4
    plot(dfana.calcVapp(CV_solutions{i}), J_values(:,2,i)) 
    hold on
end
xlim([0, 1.3])
xlabel('Voltage (V)')
ylabel('e\phi_{PL}(Acm^{-2})')
ylim([0, 2.5e-3])
legend({'10 nm', '50 nm', '90 nm', '50 nm doped'}, 'Location', 'northwest')
%% Recreate Figure 2c
figure(21)

% Work out what time value index 0.5V corresponds to
% Go from 0V to 1.3V in 120 steps => 1 step is 0.0108 mV so will get to 0.5V
% after 47 steps (to nearest integer)
time_point = 47;

% Define figure and axes
plot(CV_solutions{1}.x + 80*10^-7, log10(CV_solutions{1}.u(time_point,:,2)), 'b--')

hold on

% Do background colours
% Arguments of patch give the coordinates of the corners of the polygon to
% be shaded
patch([0 CV_solutions{3}.x(160) CV_solutions{3}.x(160) 0], [18 18 0 0], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
patch([CV_solutions{3}.x(160) CV_solutions{3}.x(260) CV_solutions{3}.x(260) CV_solutions{3}.x(160)], [18 18 0 0], 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch([CV_solutions{3}.x(660) CV_solutions{3}.x(760) CV_solutions{3}.x(760) CV_solutions{3}.x(660)], [18 18 0 0], 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch([CV_solutions{3}.x(760) CV_solutions{3}.x(end) CV_solutions{3}.x(end) CV_solutions{3}.x(760)], [18 18 0 0],  'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none')

%Plot actual data
plot(CV_solutions{1}.x + 80*10^-7, log10(CV_solutions{1}.u(time_point,:,2)), 'b--')
plot(CV_solutions{1}.x + 80*10^-7, log10(CV_solutions{1}.u(time_point,:,3)), 'r--')
plot(CV_solutions{3}.x, log10(CV_solutions{3}.u(time_point,:,2)),'b')
plot(CV_solutions{3}.x, log10(CV_solutions{3}.u(time_point,:,3)), 'r')

hold off

% Labels and tick marks
xlim([0, sum(CV_solutions{3}.par.d)])
xticklabels({'0', '100', '200', '300', '400', '500'})
ylim([8,18])
yticklabels({'10^{8}', '10^{10}', '10^{12}', '10^{14}', '10^{16}', '10^{18}'})
xlabel('distance (nm)')
ylabel('carrier concentration (cm^{-3})')
title('0.5V applied')
ax=gca;
ax.TitleHorizontalAlignment = 'left';
legend({'','','','','',' n, 10 nm', ' p, 10 nm',' n, 90 nm', ' p, 90 nm' }, 'Position', [0.7 0.1 0.2 0.2])
