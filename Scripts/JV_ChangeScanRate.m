% Make plot of JV curves at different scan rates
%prebias = [0.7, 0.8, 0.9, 1.0, stats_ion.Vocf, 1.1];
prebias = [0.55, 0.625, 0.70, 0.80, 0.90, 1.00, stats_ion.Voc_f];
scan_rates = [1e-2, 1e-1, 1, 5, 50, 100, 1000];
JV_Solutions_ScanRate = cell(length(prebias),length(scan_rates));
error_log = zeros(length(prebias),length(scan_rates));

illuminated_sol = changeLight(eqm_QJV.ion, 1.1, 0, 1);

%prebias 0.90 V didn't work
for j = 1:length(prebias)
    try
        biased_eqm = genVappStructs(illuminated_sol, prebias(j), 1);
    catch
        warning(['Could not apply prebias for ' num2str(prebias(j)) ' V.'])
    end
    for i = 1:length(scan_rates)
        try           
            disp("Reached 1 Sun bias eqm, starting JV protocol")
            JV_Solutions_ScanRate{j,i} = doCV(biased_eqm, 1.1, prebias(j), -0.20, 1.20, scan_rates(i), 1, 281);            
        catch
            warning(['JV scan failed for scan rate ' num2str(scan_rates(i)) ' Vs-1.'])
            error_log(j,i) = 1;
        end
    end
end

%% Get J for plotting 
points = length(JV_Solutions_ScanRate{1,1}.t);
Vapp = zeros(length(prebias), points);
J_values = zeros(length(prebias), length(scan_rates), points);

for j = 1:length(prebias)
    Vapp(j,:) = dfana.calcVapp(JV_Solutions_ScanRate{j,1});
    for i = 1:length(scan_rates)
        try
            Jtemp = dfana.calcJ(JV_Solutions_ScanRate{j,i});
            J_values(j,i,:) = Jtemp.tot(:,1);
        catch
            warning(['No JV solution available for scan rate of ' num2str(scan_rates(i)) ' Vs-1.'])
            J_values(j,i,:) = zeros(1,1,points);
        end
    end
end
%% Plotting 
figure('Name', 'Scan Rate Dependent JV', 'Position', [100 100 1250 1250])
num = 3;
for k = 1:length(scan_rates)
    Jplot = squeeze(J_values(num,k,:));
    plot(Vapp(num,:), 1000*Jplot)
    hold on
end

xline(0, 'Color', 'black')
yline(0, 'Color', 'black')
hold off
legend({'0.01 Vs^{-1}','0.1 Vs^{-1}', '1 Vs^{-1}', '5 Vs^{-1}', '50 Vs^{-1}', '100 Vs^{-1}','1000 Vs^{-1}'}, 'Location', 'northwest')
xlim([-0.1, 1.2])
ylim([-25, 25])
xlabel('Voltage(V)')
ylabel('Current Density (Acm^{-2})')