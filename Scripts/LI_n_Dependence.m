

%% Read in data
par1 = pc('Input_files/HTL_MAPI_NegOffset.csv');
par2 = pc('Input_files/HTL_MAPI_NegOffset.csv');
par3 = pc('Input_files/HTL_MAPI_PosOffset.csv');

pars = {par1};
%pars = {par1, par2, par3};
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

%% Store Voc Values
FF = zeros(num_devices, num_samples);
for j = 1:num_samples
    for i = 1:num_devices
    JVstats{i,j} = CVstats(JVsols{i,j});
    Voc(i,j) = JVstats{i,j}.Voc_f;
    end
end

%% Voc vs LI and FF vs LI
colours = {[0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.8500 0.3250 0.0980]};
markers = {'s', 'd', 'o'};

figure('Name', 'Voc vs LI', 'Position', [100, 100, 1000, 1000])
box on
for i = 1:num_devices
    %order = [2,1,3];
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

%% Charge Stored in different layers vs LI 
num_stop_HTL = JVsols{1,1}.par.layer_points(1);
num_start_pero = sum(JVsols{1,1}.par.layer_points(1:2))+1;
num_stop_pero = num_start_pero + JVsols{1,1}.par.layer_points(3)-1;
num_start_ETL = sum(JVsols{1,1}.par.layer_points(1:4))+1;

x = JVsols{1,1}.x;
Vapp = dfana.calcVapp(JVsols{1,1});
num_t = length(JVsols{1,1}.t);
area = 0.045; 
e = 1.6e-19;
device_num = 1;

Qn = zeros(num_devices,num_samples,3,2);
Qp = zeros(num_devices,num_samples,3,2);

for i=1:num_devices
    for j=1:num_samples
        %Dark eqm Values
        Qn(i,j,1,1) = e*area*trapz(x(1:num_stop_HTL), eqm{i}.ion.u(end,1:num_stop_HTL,2),2);
        Qp(i,j,1,1) = e*area*trapz(x(1:num_stop_HTL), eqm{i}.ion.u(end,1:num_stop_HTL,3),2);
        Qn(i,j,2,1) = e*area*trapz(x(num_start_pero:num_stop_pero), eqm{i}.ion.u(end,num_start_pero:num_stop_pero,2),2);
        Qp(i,j,2,1) = e*area*trapz(x(num_start_pero:num_stop_pero), eqm{i}.ion.u(end,num_start_pero:num_stop_pero,3),2);
        Qn(i,j,3,1) = e*area*trapz(x(num_start_ETL:end), eqm{i}.ion.u(end,num_start_ETL:end,2),2);
        Qp(i,j,3,1) = e*area*trapz(x(num_start_ETL:end), eqm{i}.ion.u(end,num_start_ETL:end,3),2);

        %Interpolate solutions at Voc
        %u has form u(time:position:species)
        stop = ceil(num_t/2);
        n = interp1(Vapp(1:stop), JVsols{i,j}.u(1:stop,:,2), Voc(j));
        p = interp1(Vapp(1:stop), JVsols{i,j}.u(1:stop,:,3), Voc(j));

        Qn(i,j,1,2) = e*area*trapz(x(1:num_stop_HTL), n(1:num_stop_HTL));
        Qp(i,j,1,2) = e*area*trapz(x(1:num_stop_HTL), p(1:num_stop_HTL));
        Qn(i,j,2,2) = e*area*trapz(x(num_start_pero:num_stop_pero), n(num_start_pero:num_stop_pero));
        Qp(i,j,2,2) = e*area*trapz(x(num_start_pero:num_stop_pero), p(num_start_pero:num_stop_pero));
        Qn(i,j,3,2) = e*area*trapz(x(num_start_ETL:end), n(num_start_ETL:end));
        Qp(i,j,3,2) = e*area*trapz(x(num_start_ETL:end), p(num_start_ETL:end));
    end
end

%%
%Plot charge vs Voc data
figure('Name', 'Q vs Voc')
box on 
semilogy(Voc, Qp(1,:,1,2), 'color', 'red', 'Marker', 'o')
hold on
semilogy(Voc, Qn(1,:,2,2), 'color', 'blue', 'Marker', 'x')
semilogy(Voc, Qp(1,:,2,2), 'color', 'red', 'Marker', 'x')
semilogy(Voc, Qn(1,:,3,2), 'color', 'blue', 'Marker', 's')
ylabel('Charge (C)')
xlabel('Open-Circuit Voltage (V)')

%%
%Just ETL electrons on linear axis
figure('Name', 'Q vs Voc')
box on 
plot(Voc, Qn(1,:,3,2), 'color', 'blue', 'Marker', 's')
ylabel('Charge (C)')
xlabel('Open-Circuit Voltage (V)')