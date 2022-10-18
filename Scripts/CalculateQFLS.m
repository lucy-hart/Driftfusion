%Calculate QFLS for Weidong ETL stuff. Run this after Weidong_ETL

%% Calculate QFLS 
num_start = sum(CV_solutions_ion{1}.par.layer_points(1:2))+1;
num_stop = num_start + CV_solutions_ion{1}.par.layer_points(3)-1;
num_values = length(CV_solutions_ion{1}.t);
x = cell(1, num_devices);
for i = 1:num_devices
    x{i} = CV_solutions_ion{i}.par.x_sub;
end
d = CV_solutions_ion{1}.par.d(3);
QFLS_ion = zeros(num_values,num_devices);
QFLS_el = zeros(num_values,num_devices);
for y=1:num_devices
    [Ecb_ion, Evb_ion, Efn_ion, Efp_ion] = dfana.calcEnergies(CV_solutions_ion{y});
    [Ecb_el, Evb_el, Efn_el, Efp_el] = dfana.calcEnergies(CV_solutions_el{y});
    QFLS_ion(:,y) = trapz(x{y}(num_start:num_stop), Efn_ion(:, num_start:num_stop)-Efp_ion(:,num_start:num_stop),2)/d;
    QFLS_el(:,y) = trapz(x{y}(num_start:num_stop), Efn_el(:,num_start:num_stop)-Efp_el(:,num_start:num_stop),2)/d;
end

QFLS_SC_ion = QFLS_ion(291,:);
QFLS_SC_el = QFLS_el(291,:);

%% Calculate OC time points
QFLS_OC_ion = zeros(1,num_devices);
QFLS_OC_el = zeros(1,num_devices);

for z=1:num_devices
    QFLS_OC_ion(z) = interp1(v(ceil(num_values/2):end), QFLS_ion(ceil(num_values/2):end,z), Voc_values_ion(z,1));
    QFLS_OC_el(z) = interp1(v(ceil(num_values/2):end), QFLS_el(ceil(num_values/2):end,z), Voc_values_el(z,1));
end

%% Find 'Figure of Merit'

Delta_mu_ion = (QFLS_OC_ion-QFLS_SC_ion)*1000;
Delta_mu_el = (QFLS_OC_el-QFLS_SC_el)*1000;
QFLS_Loss = (QFLS_OC_ion-Voc_values_ion(:,1))*1000;

%% Plot
figure('Name', 'QFLSPlot Mobility', 'Position', [100 100 1250 2000])
mob_order = [2,3,1];
for i = 1:3
    j = mob_order(i);
    plot(V(161:end), QFLS_ion(161:end,i), 'color', colors_JV{j}, 'LineWidth', 3) 
    hold on 
    %plot(Voc_values_ion(i,1), QFLS_OC_ion(i), 'ko', 'MarkerSize', 15, 'LineWidth', 2)
    plot(V(161:end), QFLS_el(161:end,i), '--', 'color', colors_JV{j}, 'LineWidth', 3) 
    %plot(Voc_values_el(i,1), QFLS_OC_el(i), 'ko', 'MarkerSize', 15, 'LineWidth', 2)
end
hold off

set(gca, 'FontSize', 25)
xlim([0, 1.3])
ylim([0.85, 1.2])
legend({'  10^{-3}', '', '  10^{-4}', '', '  10^{-5}', ''}, 'Location', 'southeast', 'FontSize', 30)
%legend({'  10^{-3}', '', '','','  10^{-4}', '', '', '', '  10^{-5}', '', '', ''}, 'Location', 'southeast', 'FontSize', 30)
xlabel('Voltage (V)', 'FontSize', 30)
ylabel('QFLS (eV)', 'FontSize', 30)
title(legend, 'ETM Mobility (cm^{2}/V/s)', 'Fontsize', 30)
ax4 = gca;

figure('Name', 'QFLSPlot Energetics', 'Position', [100 100 1250 2000])
energetics_order_devices = [4,1,5];
for j = 1:3
    i = energetics_order_devices(j);
    plot(V(161:end), QFLS_ion(161:end,i), 'color', colors_JV{j}, 'LineWidth', 3) 
    hold on 
    %plot(Voc_values_ion(i,1), QFLS_OC_ion(i), 'ko', 'MarkerSize', 15, 'LineWidth', 2)
    plot(V(161:end), QFLS_el(161:end,i), '--', 'color', colors_JV{j}, 'LineWidth', 3) 
    %plot(Voc_values_el(i,1), QFLS_OC_el(i), 'ko', 'MarkerSize', 15, 'LineWidth', 2)
end

set(gca, 'FontSize', 25)
xlim([0, 1.3])
ylim([0.85, 1.2])
xlabel('Voltage (V)', 'FontSize', 30)
ylabel('QFLS (eV)', 'FontSize', 30)
legend({'  4.15', '', '  3.95', '', '  3.75', ''}, 'Location', 'southeast', 'FontSize', 30)
%legend({'  4.15', '', '', '', '  3.95', '', '', '', '  3.75', '', '', ''}, 'Location', 'southeast', 'FontSize', 30)
title(legend, '|ETM LUMO| (eV)', 'Fontsize', 30)
ax5 = gca;

%%
figure('Name', 'QFLSPlot Energetics', 'Position', [100 100 1250 2000])
energetics_order_devices = [4,1,5];
box on
hold on 
plot(V(161:end), QFLS_ion(161:end,1), 'color', colors_JV{2}, 'LineWidth', 3) 
plot(Voc_values_ion(1,1), QFLS_OC_ion(1), 'ko', 'MarkerSize', 15, 'LineWidth', 2)
plot(V(161:end), QFLS_el(161:end,1), '--', 'color', colors_JV{2}, 'LineWidth', 3) 
plot(Voc_values_el(i,1), QFLS_OC_el(1), 'ko', 'MarkerSize', 15, 'LineWidth', 2)

set(gca, 'FontSize', 25)
xlim([0, 1.3])
ylim([0.85, 1.2])
xlabel('Voltage (V)', 'FontSize', 30)
ylabel('QFLS (eV)', 'FontSize', 30)
legend({'  Mobile Ions', '','  No Mobile Ions',''}, 'Location', 'southeast', 'FontSize', 30)
ax6 = gca;

%% Save Plots at 300 dpi
save = 1;
fig_num = 2;

if save == 1 && fig_num == 1
    exportgraphics(ax4, ...
    'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Weidong_ETL\Paper\v2\QFLS_mob_NoCircles.png', ...
    'Resolution', 300)
elseif save == 1 && fig_num == 2
    exportgraphics(ax5, ...
    'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Weidong_ETL\Paper\v2\QFLS_energetics_NoCircles.png', ...
    'Resolution', 300)
elseif save == 1 && fig_num == 3
    exportgraphics(ax6, ...
    'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Weidong_ETL\Paper\v2\QFLS.png', ...
    'Resolution', 300)
end 