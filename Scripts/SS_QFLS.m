%Script for investigating amount of QFLS at SC and quenching of PL between
%OC and SC as a function of interlayer material and width
%Parameters and CV scan rate taken from 'Voltage-Dependent
%Photoluminescence and How It Correlates with the Fill Factor 
%and Open-Circuit Voltage in Perovskite Solar Cells'.

%% Read in data and create devices with HTLs of different thicknesses
par_10 = pc('Input_files/ptaa_mapi_pcbm.csv');
par_10.vsr_mode=0;

par_10 = refresh_device(par_10);

par_50 = par_10;
par_50.d(1) = 5*10^-6;

par_50 = refresh_device(par_50);

par_90 = par_10;
par_90.d(1) = 9*10^-6;
par_90 = refresh_device(par_90);

devices = {par_10, par_50, par_90};

%% Find eqm solutions at 1 sun intensity
eqm_solutions_dark = cell(1,3);
eqm_solutions_light = cell(1,3);
for i = 1:3
    eqm = equilibrate(devices{i});
    eqm_solutions_dark{i} = eqm;
    eqm_solutions_light{i} = changeLight(eqm.ion,1,0);
end

%% Perform CV scans
CV_solutions = cell(1,3);
for j = 1:3
    sol = eqm_solutions_dark{j}.ion;
    CV_solutions{j} = doCV(sol, 1, 0, 1.3, 0, 10e-3, 1, 241);
    dfplot.JtotVapp(CV_solutions{j},0)
    hold on
end
hold off
legend({'10 nm', '50 nm', '90 nm'}, 'Location', 'northwest')
xlim([0, 1.3])
ylim([-0.03, 0.01])

%% Break down contributions to the current
%Columns in J_values are J_gen, J_rad, J_nonrad and J_ext
J_values = zeros(241,4,3);
e = -1.60*10^-19;
for k = 1:3
    w = eqm_solutions_dark{k}.el.par.d;
    num_points = eqm_solutions_dark{k}.el.par.layer_points;
    volume_matrix = zeros(241,sum(num_points));
    volume_matrix(:,1:num_points(1))= w(1)/num_points(1);
    start = num_points(1);
    for m = 2:5
        volume_matrix(:,(start+1):start+num_points(m)) = w(m)/num_points(m);
        start = start+num_points(m);
    end
    loss_currents = dfana.calcr(CV_solutions{k},'whole');
    J_values(:,1,k) = e*sum(dfana.calcg(CV_solutions{k}).*volume_matrix,2);
    J_values(:,2,k) = e*sum(loss_currents.btb(:,1:end-1).*volume_matrix,2);
    J_values(:,3,k) = e*sum(loss_currents.srh(:,1:end-1).*volume_matrix,2);
    J_values(:,4,k) = dfana.calcJ(CV_solutions{k}).tot(:,1);
    
end

%% Plot contributons to the current
figure(2)
line_colour = {[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0 0.4470 0.7410], [0.4660 0.6740 0.1880]};
V_10 = dfana.calcVapp(CV_solutions{1});
V_90 = dfana.calcVapp(CV_solutions{3});
hold on
for n = 1:4
    plot(V_10(1:121), J_values(1:121,n,3), 'color', line_colour{n}) 
    plot(V_90(1:121), J_values(1:121,n,1), 'color', line_colour{n}, 'linestyle', '--')
end
plot(V_10(1:121), zeros(1,121), 'black', 'LineWidth', 1)
hold off
xlim([0, 1.3])
xlabel('Voltage (V)')
ylim([-0.05, 0.01])
ylabel('Current Density (Acm^{-2})')
legend({'J_{gen}', '', 'J_{rad}', '', 'J_{nonrad}', '', 'J_{ext}', '', 'J = 0'}, 'Location', 'bestoutside')

%% Plot 'PL' results
figure(11)
for i = 1:3
    semilogy(dfana.calcVapp(CV_solutions{i}), -J_values(:,2,i)) 
    hold on
end
xlim([0, 1.3])
xlabel('Voltage (V)')
ylabel('e\phi_{PL}')
legend({'10 nm', '50 nm', '90 nm'}, 'Location', 'northwest')
%% Recreate Figure 2c
figure(21)
% Work out what time value index 0.5V corresponds to
% Go from 0V to 1.3V in 120 steps => 1 step is 0.0108 mV so will get to 0.5V
% after 47 steps (to nearest integer)

% Define figure and axes
plot(CV_solutions{1}.x + 80*10^-7, log10(CV_solutions{1}.u(47,:,2)), 'b--')

hold on

% Do background colours
% Arguments of patch give the coordinates of the corners of the polygon to
% be shaded
patch([0 CV_solutions{3}.x(160) CV_solutions{3}.x(160) 0], [18 18 0 0], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
patch([CV_solutions{3}.x(160) CV_solutions{3}.x(260) CV_solutions{3}.x(260) CV_solutions{3}.x(160)], [18 18 0 0], 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch([CV_solutions{3}.x(660) CV_solutions{3}.x(760) CV_solutions{3}.x(760) CV_solutions{3}.x(660)], [18 18 0 0], 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch([CV_solutions{3}.x(760) CV_solutions{3}.x(end) CV_solutions{3}.x(end) CV_solutions{3}.x(760)], [18 18 0 0],  'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none')

%Plot actual data
plot(CV_solutions{1}.x + 80*10^-7, log10(CV_solutions{1}.u(47,:,2)), 'b--')
plot(CV_solutions{1}.x + 80*10^-7, log10(CV_solutions{1}.u(47,:,3)), 'r--')
plot(CV_solutions{3}.x, log10(CV_solutions{3}.u(47,:,2)),'b')
plot(CV_solutions{3}.x, log10(CV_solutions{3}.u(47,:,3)), 'r')

hold off

% Labels and tick marks
xlim([0, sum(CV_solutions{3}.par.d)])
xticklabels({'0', '100', '200', '300', '400', '500'})
ylim([8,18])
yticklabels({'10^{8}', '10^{10}', '10^{12}', '10^{14}', '10^{16}', '10^{18}'})
xlabel('distance (nm)')
ylabel('carrier concentration (cm^{-3})')
legend({'','','','','',' n, 10 nm', ' p, 10 nm',' n, 90 nm', ' p, 90 nm' }, 'Position', [0.7 0.1 0.2 0.2])
