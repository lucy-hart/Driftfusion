% Make plot of JV curves at different scan rates
%prebias = [0.7, 0.8, 0.9, 1.0, stats_ion.Vocf, 1.1];
prebias = [0.6, stats_ion.Voc_f];
scan_rates = [1e-2, 1e-1, 1, 10, 100];
JV_Solutions_ScanRate = cell(length(prebias),length(scan_rates));

Vmax = 1.2;
Vmin = -0.2;
Vdiff = Vmax - Vmin;
tsweep = Vdiff./scan_rates;

j=1;

for i = 1e-2%scan_rates
    for k = 1%:length(prebias)
        try
            illuminated_sol = changeLight(eqm_QJV.ion, 1.1, 120, 1);
            biased_eqm = VappFunction(illuminated_sol, 'sweepAndStill', [0, prebias(k), 10], 200, 100, 0);
            disp("Reached 1 Sun bias eqm, starting JV protocol")
            JV_Solutions_ScanRate{k,j} = VappFunction(biased_eqm, 'tri', [prebias(k), Vmin, Vmax, 1, tsweep(i)], tsweep(i), 281, 0);
        catch
            warning(['JV scan failed for scan rate ' num2str(i) ' Vs-1.'])
        end
        j = j+1;
    end
end

%% Plotting 
figure('Name', 'Scan Rate Dependent JV', 'Position', [100 100 1250 1250])
num = 1;
for k = 1:length(prebias)
    J = dfana.calcJ(JV_Solutions_ScanRate_fw{num,k});
    plot(v(1:151), 1000*J.tot(1:151,1))
    hold on
    plot(v(151:end), 1000*J.tot(151:end,1), 'LineStyle', '--')
    hold on
end

plot(v(:), zeros(1,length(v)), 'black', 'LineWidth', 1)
hold off
legend({'0.01 Vs^{-1}','','0.1 Vs^{-1}','', '1 Vs^{-1}','', '10 Vs^{-1}','', '100 Vs^{-1}',''}, 'Location', 'northwest')
xlim([-0.1, 1.2])
ylim([-25, 5])
xlabel('Voltage(V)')
ylabel('Current Density (Acm^{-2})')