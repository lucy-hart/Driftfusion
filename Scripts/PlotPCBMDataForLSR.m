data = readtable("QFLS_vs_Vap_PCBM_Experimental.xlsx", 'Range', 'A4:B134');

%%
figure('Name', 'QFLSPlotExperimental', 'Position', [100 100 1250 2000])
box on
hold on 
plot(data{:,1}, data{:,2}, 'color', [0.4660 0.6740 0.1880], 'LineWidth', 3) 


set(gca, 'FontSize', 35)
xlim([0, 1.3])
ylim([1, 1.2])
xlabel('Voltage (V)', 'FontSize', 40)
ylabel('QFLS (eV)', 'FontSize', 40)
%legend({'  Mobile Ions', '','  No Mobile Ions',''}, 'Location', 'southeast', 'FontSize', 40)
ax1 = gca;

%%
save = 1;

if save == 1 
    exportgraphics(ax1, ...
    'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\LSR\QFLS_PCBM_Experimental.png', ...
    'Resolution', 300)
end