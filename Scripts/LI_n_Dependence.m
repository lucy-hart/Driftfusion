% See how distribution of carriers in PCBM and ICBA devices varies with intensity

%% Read in data
par1 = pc('Input_files/PTAA_MAPI_PCBM_v5.csv');
par2 = pc('Input_files/PTAA_MAPI_PCBM_HigherLUMO.csv');
par3 = pc('Input_files/PTAA_MAPI_PCBM_LowerLUMO.csv');

pars = {par1, par2, par3};
num_devices = length(pars);
for i = 1:num_devices
    par = pars{i};
    par.light_source1 = 'laser';
    par.laser_lambda1 = 532;
    par.pulsepow = 62;
    par.RelTol_vsr = 0.1;
    pars{i} = refresh_device(par);
end

eqm = cell(1,num_devices);

for i = 1:num_devices
    eqm{i} = equilibrate(pars{i});
end
%% Set up parameters 
num_samples = 10;

LightInt = logspace(log10(0.1), log10(5), num_samples);
Voc = zeros(num_devices, num_samples);
QFLS_SC = zeros(num_devices, num_samples);
QFLS_OC = zeros(num_devices, num_samples);
JVsols = cell(num_devices, num_samples);
JVstats = cell(num_devices, num_samples);

%% Do JV Scans
error = zeros(num_devices, num_samples);
for j = 1:num_samples
    for i = 1:num_devices
        if LightInt(j) > 1.25
            try
                JVsols{i,j} = doCV(eqm{i}.ion, LightInt(j), -0.3, 1.2, -0.3, 10e-3, 1, 301);            
            catch
                warning('Could not complete JV scan. Assigning value of 0')
                JVsols{i,j} = 0;
                error(i,j) = 1;
            end
        elseif LightInt(j) <= 1.25 && LightInt(j) > 0.75
            try
                JVsols{i,j} = doCV(eqm{i}.ion, LightInt(j), -0.3, 1.1, -0.3, 10e-3, 1, 281);
            catch
                warning('Could not complete JV scan. Assigning value of 0')
                JVsols{i,j} = 0;
                error(i,j) = 1;
            end
        elseif LightInt(j) <= 0.75
            try
                JVsols{i,j} = doCV(eqm{i}.ion, LightInt(j), -0.3, 1.05, -0.3, 10e-3, 1, 271);
            catch
                warning('Could not complete JV scan. Assigning value of 0')
                JVsols{i,j} = 0;
                error(i,j) = 1; 
            end
        end
    end
end

%% Store Voc and QFLS values
num_start = sum(JVsols{1,1}.par.layer_points(1:2))+1;
num_stop = num_start + JVsols{1,1}.par.layer_points(3)-1;
x = JVsols{1,1}.par.x_sub;
d = JVsols{1,1}.par.d(3);
for j = 1:num_samples
    for i = 1:num_devices
        JVstats{i,j} = CVstats(JVsols{i,j});
        Voc(i,j) = JVstats{i,j}.Voc_r;
        num_values = length(JVsols{i,j}.t);
        [~, ~, Efn_ion, Efp_ion] = dfana.calcEnergies(JVsols{i,j});
        QFLS_ion = trapz(x(num_start:num_stop), Efn_ion(:, num_start:num_stop)-Efp_ion(:,num_start:num_stop),2)/d;
        Vapp = dfana.calcVapp(JVsols{i,j});
        start_point = ceil(num_values/2);
        QFLS_SC(i,j) = interp1(Vapp(start_point:end), QFLS_ion(start_point:end), 0);
        QFLS_OC(i,j) = interp1(Vapp(start_point:end), QFLS_ion(start_point:end), Voc(i,j));
    end
end

%% Get gradients
alpha = zeros(3, num_samples, num_devices);
kT = 0.0257;
for k = 1:num_devices
    alpha(1,:,k) = (1/kT)*gradient(Voc(k,:), log(LightInt(2))-log(LightInt(1)));
    alpha(2,:,k) = (1/kT)*gradient(QFLS_OC(k,:), log(LightInt(2))-log(LightInt(1)));
    alpha(3,:,k) = (1/kT)*gradient(QFLS_SC(k,:), log(LightInt(2))-log(LightInt(1)));
end

%% Voc vs LI and QFLS_SC/OC vs LI
colours = {[0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.8500 0.3250 0.0980]};
markers = {'s', 'd', 'o'};

figure('Name', 'Voc vs LI', 'Position', [100, 100, 1000, 1000])
box on
for i = 1:num_devices
    semilogx(LightInt, Voc(i,:), 'color', colours{i}, 'marker', markers{i}, 'MarkerFaceColor', colours{i}, ...
        'MarkerSize', 10, 'LineWidth', 2)
    if i == 1
        hold on
    end
end

xlim([0.09,6])
ylim([0.9, 1.15])
ax = gca;
ax.FontSize = 25;
xlabel('Light Intensity (% of 1 Sun)', 'FontSize', 30)
ylabel('V_{OC} (V)', 'FontSize', 30)
legend({'  3.95', '  3.75', '  4.15'}, 'Location', 'southeast', 'FontSize', 30)
title(legend, '|ETM LUMO| (eV)', 'Fontsize', 30)

figure('Name', 'QFLS vs LI', 'Position', [100, 100, 1000, 1000])
box on
for i = 1:num_devices
    semilogx(LightInt, QFLS_SC(i,:), 'color', colours{i}, 'marker', markers{i}, ...
        'MarkerSize', 10, 'LineWidth', 2)
    if i == 1
        hold on
    end
    semilogx(LightInt, QFLS_OC(i,:), 'color', colours{i}, 'marker', markers{i}, 'MarkerFaceColor', colours{i}, ...
        'MarkerSize', 10, 'LineWidth', 2)
end

xlim([0.09,6])
ylim([0.8, 1.2])
ax = gca;
ax.FontSize = 25;
xlabel('Light Intensity (% of 1 Sun)', 'FontSize', 30)
ylabel('QFLS (eV)', 'FontSize', 30)
legend({'  3.95', '', '  3.75', '', '  4.15', ''}, 'Location', 'southeast', 'FontSize', 30)
title(legend, '|ETM LUMO| (eV)', 'Fontsize', 30)

%% Get 1 Sun solutions
OneSunJVs = cell(1,num_devices);
OneSunJV_Stats = cell(1,num_devices);
for i = 1:num_devices
    OneSunJVs{i} = doCV(eqm{i}.ion, 1, -0.3, 1.2, -0.3, 10e-3, 1, 301); 
    OneSunJV_Stats{i} = CVstats(OneSunJVs{i});
end 

%% JVs at 1 Sun and 0.1 Sun
figure('Name', 'JV 1 Sun')
box on
hold on
for i = 1:num_devices
    J = dfana.calcJ(OneSunJVs{i}).tot;
    V = dfana.calcVapp(OneSunJVs{i});
    num_points = length(V);
    plot(V(1:ceil(num_points/2)), 1e3*J(1:ceil(num_points/2),1), 'color', colours{i}, 'LineWidth', 2)
    J = dfana.calcJ(JVsols{i,1}).tot;
    V = dfana.calcVapp(JVsols{i,1});
    num_points = length(V);
    plot(V(1:ceil(num_points/2)), 1000*J(1:ceil(num_points/2),1), 'color', colours{i}, 'LineWidth', 2)
end
xline(0)
yline(0)
hold off 

xlim([-0.2, 1.2])
ylim([-25, 1])
xlabel('Voltage (V)')
ylabel('Current Density (mAcm^{-2})')
legend('PCBM', '', 'ICBA', '')

