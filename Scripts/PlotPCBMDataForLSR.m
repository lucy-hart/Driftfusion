data = readtable("QFLS_vs_Vap_Experimental.xlsx", 'Range', 'A4:E135');

%%
figure('Name', 'QFLSPlotExperimental', 'Position', [100 100 1250 2000])
colours = {[0.4660 0.6740 0.1880],[0 0.4470 0.7410],[0.9290 0.6940 0.1250],[0.8500 0.3250 0.0980],};
box on
hold on 
for i = 1:4
plot(data{:,1}, data{:,i+1}, 'color', colours{i}, 'LineWidth', 3) 
end


set(gca, 'FontSize', 35)
xlim([0, 1.3])
ylim([1.03, 1.2])
xlabel('Voltage (V)', 'FontSize', 40)
ylabel('QFLS (eV)', 'FontSize', 40)
legend({'  PCBM', '  ICBA','  IPH','  KLOC-6'}, 'Location', 'northwest', 'FontSize', 40)
title(legend, 'ETL Material', 'Fontsize', 40)
ax1 = gca;

%%
save = 1;

if save == 1 
    exportgraphics(ax1, ...
    'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\LSR\QFLS_Experimental.png', ...
    'Resolution', 300)
end