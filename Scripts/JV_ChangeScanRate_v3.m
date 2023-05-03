%Make plot Voc vs prebias for fastest JV scan
%NB: give prebias values to 10 mV or Vapp won't be sampled every 10 mV in
%the JV scans
%Use QuickJVChecker.m to find Voc/flat ion potential

%% Read in data and set up parameters
par1 = pc('Input_files/PTAA_MAPI_NegOffset.csv');
par2 = pc('Input_files/PTAA_MAPI_NoOffset.csv');

devices = {par1, par2};
num_devices = length(devices);

dark_eqm = cell(1,num_devices);

for i = 1:num_devices
    dark_eqm{i} = equilibrate(devices{i});
end

prebias = [[0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15];
            [0.49, 0.60, 0.71, 0.80, 0.90, 1.00, 1.10, 1.15]];

scan_rates = [1e-2, 1000];

JV_Solutions_ScanRate_PartOne = cell(length(prebias),length(scan_rates),num_devices);
JV_Solutions_ScanRate_PartTwo = cell(length(prebias),length(scan_rates),num_devices);
JV_Solutions_ScanRate = cell(length(prebias),length(scan_rates),num_devices);
error_log = zeros(length(prebias),length(scan_rates),num_devices, 3);
num_tpoints = zeros(length(prebias),2,num_devices);

Vstart = -0.20;
Vend = 1.20;

for k = 1:num_devices
    for j = 1:length(prebias)
        num_tpoints(j,1,k) = cast(200*(prebias(k,j) - Vstart)+1, 'int32');
        num_tpoints(j,2,k) = cast(200*(Vend - prebias(k,j))+1, 'int32');
    end
end
%% Do scan rate dependent JVs

for k = 1:num_devices
    illuminated_sol = changeLight(dark_eqm{k}.ion, 1.1, 0, 1);
    for j = 1:length(prebias)
        try
            biased_eqm = genVappStructs(illuminated_sol, prebias(k,j), 1);
        catch
            warning(['Could not apply prebias for ' num2str(prebias(k,j)) ' V.'])
            continue
        end
        for i = 1:length(scan_rates)
            try           
                disp("Reached 1 Sun bias eqm, starting JV protocol")
                n = 1;
                JV_Solutions_ScanRate_PartOne{j,i,k} = doCV(biased_eqm, 1.1, prebias(k,j), Vstart, prebias(k,j), scan_rates(i), 1, num_tpoints(j,1,k));
                n = 2;
                JV_Solutions_ScanRate_PartTwo{j,i,k} = doCV(biased_eqm, 1.1, prebias(k,j), Vend, prebias(k,j), scan_rates(i), 1, num_tpoints(j,2,k));    
                n = 3;
                JV_Solutions_ScanRate{j,i,k} = doCV(biased_eqm, 1.1, prebias(k,j), Vend, Vstart, scan_rates(i), 1, 281);            
            catch
                warning(['JV scan ' num2str(n) ' failed for scan rate ' num2str(scan_rates(i)) ' Vs-1.'])
                error_log(j,i,k,n) = 1;
            end
        end
    end
end

%% Get J for plotting 
%cont for continuous and discont for discontinuous
Vapp_cont = cell(num_devices, length(prebias));
J_values_cont = cell(num_devices, length(prebias), length(scan_rates));

Vapp = cell(num_devices, length(prebias));
J_values = cell(num_devices, length(prebias), length(scan_rates));

Voc_fast = zeros(num_devices,length(prebias));

%This is convoluted but it's needed to make the JV scans continuous when you do
%a fast scan starting at a prebias voltage below Voc
for k = 1:num_devices
    for j = 1:length(prebias)
        stop1 = ceil(num_tpoints(j,1,k)/2);
        stop2 = ceil(num_tpoints(j,2,k)/2);
        V1 = dfana.calcVapp(JV_Solutions_ScanRate_PartOne{j,1,k});
        V2 = dfana.calcVapp(JV_Solutions_ScanRate_PartTwo{j,1,k});
        V = dfana.calcVapp(JV_Solutions_ScanRate{j,1,k});
        Vapp_cont{k,j} = [flip(V1(1:stop1)), V2(2:stop2)];
        start = cast(100*(Vend - prebias(k,j))+1, 'int32');
        Vapp{k,j} = V(start:start+140);
        for i = 1:length(scan_rates)
            try
                Jtemp1 = dfana.calcJ(JV_Solutions_ScanRate_PartOne{j,i,k});
                Jtemp2 = dfana.calcJ(JV_Solutions_ScanRate_PartTwo{j,i,k});
                J_values_cont{k,j,i} = [flip(Jtemp1.tot(1:stop1,1)); Jtemp2.tot(2:stop2,1)];
                Voc_fast(k,j) = interp1(J_values_cont{k,j,i}, Vapp_cont{k,j}, 0);
            catch
                warning(['No JV solution available for scan rate of ' num2str(scan_rates(i)) ' Vs-1.'])
                J_values_cont{k,j,i} = 0;
            end
            try
                Jtemp = dfana.calcJ(JV_Solutions_ScanRate{j,i,k});
                J_values{k,j,i} = Jtemp.tot(start:start+140,1); 
            catch
                warning(['No JV solution available for scan rate of ' num2str(scan_rates(i)) ' Vs-1.'])
                J_values{k,j,i} = 0;
            end
        end
    end
end
%% Plotting JV scans
%choose prebias
num = 7;
num_device = 2;

figure('Name', 'Scan Rate Dependent JV', 'Position', [50 50 800 900])
Colours = {[0 0.4470 0.7410], [0.4660 0.6740 0.1880], [0.8500 0.3250 0.0980]};
plot_fw = 1;
zero = zeros(50,1);
if plot_fw == 0
    labels = {'','',' 0.01 Vs^{-1}',' 1000 Vs^{-1}'};
elseif plot_fw == 1
    labels = {'','',' 0.01 Vs^{-1}','',' 1000 Vs^{-1}',''};
end

box on
hold on 
plot(zero(:,1), linspace(-25,25,50), 'Color', 'black', 'LineWidth', 1.5)
plot(linspace(-0.2,1.2,50), zero(:,1), 'Color', 'black', 'LineWidth', 1.5)
for k = 1:length(scan_rates)
    plot(Vapp_cont{num_device,num}, 1000*J_values_cont{num_device,num,k}, 'Color', Colours{k}, 'LineWidth', 2)
    if plot_fw == 1
        plot(Vapp{num_device,num}, 1000*J_values{num_device,num,k}, 'LineStyle', '--', 'Color', Colours{k}, 'LineWidth', 2)
    end
end

hold off
set(gca, 'Fontsize', 25)
legend(labels, 'Location', 'northwest', 'FontSize', 23)
xlim([-0.1, 1.2])
ylim([-25, 25])
xlabel('Voltage(V)', 'FontSize', 30)
ylabel('Current Density (Acm^{-2})', 'FontSize', 30)

%% Plotting Voc values versus prebias 
%Voc value of fast scan as the value in the slow scan is the same,
%irrespective of prebais
figure('Name', 'Voc versus Prebias', 'Position', [50 50 800 900])

box on 
hold on
for k = 2%1:num_devices
    plot(prebias(k,:), Voc_fast(k,:), 'Marker', 'x', 'MarkerSize', 15, 'LineStyle', 'none')
end
hold off

set(gca, 'Fontsize', 25)
xlim([0,1.2])
ylim([0,1.2])
xlabel('Prebias Voltage (V)', 'FontSize', 30)
ylabel('Open Circuit Voltage (V)', 'FontSize', 30)
%legend({'ETM 1', 'ETM 2'}, 'FontSize', 25, 'Location', 'southeast')