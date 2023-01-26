% Make plot of JV curves at different scan rates
%NB: give prebias values to 10 mV or Vapp won't be sampled every 10 mV in
%the JV scans
%Use QuickJVChecker.m to get the eqm solution and to find Voc/flat ion potential

%% Set up variables and find dark and one sun eqm solutions
%par = pc('Input_files/PTAA_MAPI_NegOffset.csv');
par = pc('Input_files/PTAA_MAPI_NoOffset.csv');

dark_eqm = equilibrate(par);

%For NegOffset case
%Flat ion potential around 0.62-0.63 meV
%Voc is 1.06 V
%prebias = [0.55, 0.63, 0.75, 0.85, 0.95, 1.06];

%For NoOffset Case
%Flat ion potential around 0.90 meV
%Voc is 1.09 V
prebias = [0.75, 0.85, 0.90, 0.95, 1.09];

scan_rates = [1e-2, 1e-1, 2e-1, 1, 2, 5, 100, 1000];

JV_Solutions_ScanRate_PartOne = cell(length(prebias),length(scan_rates));
JV_Solutions_ScanRate_PartTwo = cell(length(prebias),length(scan_rates));
JV_Solutions_ScanRate = cell(length(prebias),length(scan_rates));
error_log = zeros(length(prebias),length(scan_rates),3);
num_tpoints = zeros(length(prebias),2);

Vstart = -0.20;
Vend = 1.20;

for j = 1:length(prebias)
    num_tpoints(j,1) = cast(200*(prebias(j) - Vstart)+1, 'int32');
    num_tpoints(j,2) = cast(200*(Vend - prebias(j))+1, 'int32');
end

%One SUn slightly higher than one Sun as n and k values are a bit off
illuminated_sol = changeLight(dark_eqm.ion, 1.1, 0, 1);

%% Do the JV scans as a function of scan rate and prebias 

for j = 1:length(prebias)
    try
        biased_eqm = genVappStructs(illuminated_sol, prebias(j), 1);
    catch
        warning(['Could not apply prebias for ' num2str(prebias(j)) ' V.'])
        continue
    end
    for i = 1:length(scan_rates)
        try           
            disp("Reached 1 Sun bias eqm, starting JV protocol")
            n = 1;
            JV_Solutions_ScanRate_PartOne{j,i} = doCV(biased_eqm, 1.1, prebias(j), Vstart, prebias(j), scan_rates(i), 1, num_tpoints(j,1));
            n = 2;
            JV_Solutions_ScanRate_PartTwo{j,i} = doCV(biased_eqm, 1.1, prebias(j), Vend, prebias(j), scan_rates(i), 1, num_tpoints(j,2));    
            n = 3;
            JV_Solutions_ScanRate{j,i} = doCV(biased_eqm, 1.1, prebias(j), Vstart, Vend, scan_rates(i), 1, 281);            
        catch
            warning(['JV scan ' num2str(n) ' failed for scan rate ' num2str(scan_rates(i)) ' Vs-1.'])
            error_log(j,i,n) = 1;
        end
    end
end

%% Get J for plotting 
%cont for continuous and discont for discontinuous
Vapp_cont = cell(1,length(prebias));
J_values_cont = cell(length(prebias), length(scan_rates));

Vapp = cell(1,length(prebias));
J_values = cell(length(prebias), length(scan_rates));

Voc_fast = zeros(1,length(prebias));

%This is convoluted but it's needed to make the JV scans continuous when you do
%a fast scan starting at a prebias voltage below Voc
for j = 1:length(prebias)
    stop1 = ceil(num_tpoints(j,1)/2);
    stop2 = ceil(num_tpoints(j,2)/2);
    V1 = dfana.calcVapp(JV_Solutions_ScanRate_PartOne{j,1});
    V2 = dfana.calcVapp(JV_Solutions_ScanRate_PartTwo{j,1});
    V = dfana.calcVapp(JV_Solutions_ScanRate{j,1});
    Vapp_cont{1,j} = [flip(V1(1:stop1)), V2(2:stop2)];
    start = cast(100*(prebias(j)-Vstart)+1, 'int32');
    Vapp{1,j} = V(start:start+140);
    for i = 1:length(scan_rates)
        try
            Jtemp1 = dfana.calcJ(JV_Solutions_ScanRate_PartOne{j,i});
            Jtemp2 = dfana.calcJ(JV_Solutions_ScanRate_PartTwo{j,i});
            J_values_cont{j,i} = [flip(Jtemp1.tot(1:stop1,1)); Jtemp2.tot(2:stop2,1)];
            Voc_fast(j) = interp1(J_values_cont{j,i}, Vapp{1,j}, 0);
        catch
            warning(['No JV solution available for scan rate of ' num2str(scan_rates(i)) ' Vs-1.'])
            J_values_cont{j,i} = 0;
        end
        try
            Jtemp = dfana.calcJ(JV_Solutions_ScanRate{j,i});
            J_values{j,i} = Jtemp.tot(start:start+140,1); 
        catch
            warning(['No JV solution available for scan rate of ' num2str(scan_rates(i)) ' Vs-1.'])
            J_values{j,i} = 0;
        end
    end
end
%% Plotting 
figure('Name', 'Scan Rate Dependent JV', 'Position', [50 50 800 900])
%choose prebias
num =  3;
Colours = parula(length(scan_rates));
plot_fw = 1;
plot_bw = 1;
zero = zeros(50,1);
if plot_fw == 0
    labels = {'','',' 0.01 Vs^{-1}',' 0.1 Vs^{-1}',' 0.2 Vs^{-1}',' 1 Vs^{-1}',' 2 Vs^{-1}',' 5 Vs^{-1}',' 100 Vs^{-1}',' 1000 Vs^{-1}'};
elseif plot_fw == 1
    labels = {'','',' 0.01 Vs^{-1}','',' 0.1 Vs^{-1}','', ' 0.2 Vs^{-1}', '',' 1 Vs^{-1}', '',' 2 Vs^{-1}', '',' 5 Vs^{-1}', '',' 100 Vs^{-1}','',' 1000 Vs^{-1}',''};
end

box on
hold on 
plot(zero(:,1), linspace(-25,25,50), 'Color', 'black', 'LineWidth', 1.5)
plot(linspace(-0.2,1.2,50), zero(:,1), 'Color', 'black', 'LineWidth', 1.5)
for k = 1:length(scan_rates)
    if plot_bw == 1
        plot(Vapp_cont{num}, 1000*J_values_cont{num,k}, 'Color', Colours(k,:), 'LineWidth', 2)
    end
    if plot_fw == 1
        plot(Vapp{num}, 1000*J_values{num,k}, 'LineStyle', '--', 'Color', Colours(k,:), 'LineWidth', 2)
    end
end

hold off
set(gca, 'Fontsize', 25)
legend(labels, 'Location', 'northwest', 'FontSize', 23)
xlim([-0.1, 1.2])
ylim([-25, 25])
xlabel('Voltage(V)', 'FontSize', 30)
ylabel('Current Density (Acm^{-2})', 'FontSize', 30)