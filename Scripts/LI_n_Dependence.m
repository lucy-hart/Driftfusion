% See how distribution of carriers in PCBM and ICBA devices varies with intensity

%% Read in data
par1 = pc('Input_files/HTL_MAPI_NoOffset.csv');
par2 = pc('Input_files/HTL_MAPI_NegOffset.csv');
par3 = pc('Input_files/HTL_MAPI_PosOffset.csv');

pars = {par1, par2, par3};
num_devices = length(pars);
eqm = cell(1,num_devices);

for i = 1:num_devices
    eqm{i} = equilibrate(pars{i});
end
%% Set up parameters 
num_samples = 10;

LightInt = logspace(log10(0.1), log10(5), num_samples);
Voc = zeros(num_devices, num_samples);
OC_time_idx = zeros(num_devices, num_samples);
JVsols = cell(num_devices, num_samples);
JVstats = cell(num_devices, num_samples);

%% Do JV Scans
error = zeros(num_devices, num_samples);
for j = 1:num_samples
    for i = 1:num_devices
        if LightInt(j) > 1.5
            try
                JVsols{i,j} = doCV(eqm{i}.ion, LightInt(j), -0.3, 1.4, -0.3, 10e-3, 1, 341);            
            catch
                warning('Could not complete JV scan. Assigning value of 0')
                JVsols{i,j} = 0;
                error(i,j) = 1;
            end
        elseif LightInt(j) <= 1.5 && LightInt(j) > 0.75
            try
                JVsols{i,j} = doCV(eqm{i}.ion, LightInt(j), -0.3, 1.3, -0.3, 10e-3, 1, 321);
            catch
                warning('Could not complete JV scan. Assigning value of 0')
                JVsols{i,j} = 0;
                error(i,j) = 1;
            end
        elseif LightInt(j) <= 0.75
            try
                JVsols{i,j} = doCV(eqm{i}.ion, LightInt(j), -0.3, 1.2, -0.3, 10e-3, 1, 301);
            catch
                warning('Could not complete JV scan. Assigning value of 0')
                JVsols{i,j} = 0;
                error(i,j) = 1; 
            end
        end
    end
end

%% Store Voc and FF Values Values
FF = zeros(num_devices, num_samples);
for j = 1:num_samples
    for i = 1:num_devices
    JVstats{i,j} = CVstats(JVsols{i,j});
    Voc(i,j) = JVstats{i,j}.Voc_f;
    %FF(i,j) = JVstats{i,j}.FF_f;
    end
end

%% Voc vs LI and FF vs LI
colours = {[0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.8500 0.3250 0.0980]};
markers = {'s', 'd', 'o'};

figure('Name', 'Voc vs LI', 'Position', [100, 100, 1000, 1000])
box on
for i = 1:num_devices
    semilogx(LightInt, (1/0.0257)*gradient(Voc(i,:), log(LightInt(:))) , 'color', colours{i}, 'marker', markers{i}, 'MarkerFaceColor', colours{i}, ...
        'MarkerSize', 10, 'LineWidth', 2)
    if i == 1
        hold on
    end
end

xlim([0.09,6])
ylim([1, 2])
ax = gca;
ax.FontSize = 25;
xlabel('Light Intensity (% of 1 Sun)', 'FontSize', 30)
ylabel('Suns-V_{OC} Ideality', 'FontSize', 30)
title(legend, 'E_{LUMO} - E_{CB}')
legend({' -0.2 eV', '  0.0 eV', '+0.2 eV'}, 'Location', 'northwest', 'FontSize', 30)
%%
figure('Name', 'FF vs LI', 'Position', [100, 100, 1000, 1000])
box on
for i = 1:num_devices
    semilogx(LightInt, FF(i,:), 'color', colours{i}, 'marker', markers{i}, 'MarkerFaceColor', colours{i}, ...
        'MarkerSize', 10,'LineWidth', 2)
    if i == 1
        hold on
    end
end

xlim([0.09,6])
ylim([0.65, 0.8])
ax = gca;
ax.FontSize = 25;
xlabel('Light Intensity (% of 1 Sun)', 'FontSize', 30)
ylabel('Fill Factor', 'FontSize', 30)
legend('PCBM', 'ICBA', 'Location', 'southwest', 'FontSize', 30)

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

