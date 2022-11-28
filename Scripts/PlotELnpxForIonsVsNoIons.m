%Plot parameters
lw = 3;
fontsize = 30;
ticksize = 25;
save = 1;

%% Find OC times
num_devices = 3;
Voc_ion = zeros(1,num_devices);
Voc_el = zeros(1,num_devices);
OC_time_ion = zeros(1,num_devices);
OC_time_el = zeros(1,num_devices);

for z=1:num_devices
    V_temp = dfana.calcVapp(JV_solutions_ion{z});
    Voc_ion(z) = CVstats(JV_solutions_ion{z}).Voc_f;
    OC_time_ion(z) = find(abs(Voc_ion(z)-V_temp) == min(abs(Voc_ion(z)-V_temp)),1);
    Voc_el(z) = CVstats(JV_solutions_el{z}).Voc_f;
    OC_time_el(z) = find(abs(Voc_el(z)-V_temp) == min(abs(Voc_el(z)-V_temp)),1);
end

%% Plot E field at SC with and without effects of mobile ions
%and ion density at SC
figure('Name', 'Field Screening', 'Position', [100 100 1250 2000])
set(gca, 'DefaultLineLineWidth', lw, 'FontSize', ticksize)

num = 2;
if num ~= 1
    pos_ion = [0.66 0.8 0.1 0.05];
else
    pos_ion = [0.68 0.8 0.1 0.05];
end

upper_lim = [3.5e4, 4.5e4, 4.5e4];
ytick_cell = {[0 0.5e4 1.0e4 1.5e4 2.0e4 2.5e4, 3.0e4];
    [0  1.0e4  2.0e4  3.0e4  4.0e4];
    [0 1.0e4 2.0e4 3.0e4 4.0e4]};

T = 31;

FV_ion = dfana.calcF(JV_solutions_ion{num}, 'whole');
FV_el = dfana.calcF(JV_solutions_el{num}, 'whole');
x_values = JV_solutions_ion{num}.x * 1e7;

hold on

plot(x_values-12, -FV_ion(T,:), 'Color', [0.4940 0.1840 0.5560])
plot(x_values-12, -FV_el(T,:), 'Color', [0.4940 0.1840 0.5560], 'LineStyle', '--')
plot(x_values-12, zeros(length(x_values)), 'Color', 'black', 'LineWidth', 2)

hold off

box on
xlim([0, 400])
ylim([-0.5e4, upper_lim(num)])
xticks([0 100 200 300 400])
yticks(ytick_cell{num})
ylabel('|Electric Field Strength| (Vcm^{-1})', 'FontSize', fontsize)
xlabel('Perovskite depth (nm)', 'FontSize', fontsize)
legend('F_{ion}', 'F_{no ions}', '', 'NumColumns', 1, 'Position', [0.72 0.8 0.15 0.02], 'FontSize', fontsize)

fig1a = gcf;

%Ions SC
upper_lim_ion = [1.4e18, 1.5e18, 1.4e18];

figure('Name', 'Cation Densities SC', 'Position', [100 100 1250 2000])
set(gca, 'DefaultLineLineWidth', lw, 'FontSize', ticksize)

hold on

% Do background colours
% Arguments of patch give the coordinates of the corners of the polygon to
% be shaded
patch('XData',[0 x_values(500) x_values(500) 0], 'YData',[1e19 1e19 1e9 1e9], 'FaceColor', colors_JV{1},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
patch('XData',[x_values(500) x_values(600) x_values(600) x_values(500)], 'YData',[1e19 1e19 1e9 1e9], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1100) x_values(1200) x_values(1200) x_values(1100)], 'YData',[1e19 1e19 1e9 1e9], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1200) x_values(end) x_values(end) x_values(1200)], 'YData',[1e19 1e19 1e9 1e9],  'FaceColor', colors_JV{3},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')

T = 31;

semilogy(x_values, (JV_solutions_ion{num}.u(T,:,4)), 'color', colors_JV{2}, 'LineStyle', '-')

hold off 

ylim([1e8, upper_lim_ion(num)])
if num == 1
    xlim([0, x_values(1540)])
else
    xlim([0, max(x_values)])
end
ylabel('Cation Density (cm^{-3})', 'Fontsize', fontsize)
xlabel('Device depth (nm)', 'Fontsize', fontsize)
xticks([0 100 200 300 400])
if num == 2
    yticks([2e17 4e17 6e17 8e17 10e17 12e17 14e17])
else
    yticks([2e17 4e17 6e17 8e17 10e17 12e17])
end
legend('','','','','Short Circuit', 'NumColumns', 1, 'Position', pos_ion, 'Fontsize', fontsize)
title(legend, 'Ion Distribution', 'Fontsize', fontsize)

fig1b = gcf;

%% Plot E field at OC with and without effects of mobile ions
figure('Name', 'Field at OC', 'Position', [100 100 1250 2000])
set(gca, 'DefaultLineLineWidth', lw, 'FontSize', ticksize)

lower_lim = [3.5e4, -3e5, 4.5e4];
upper_lim = [3.5e4, 3e5, 4.5e4];
ytick_cell = {[0 0.5e4 1.0e4 1.5e4 2.0e4 2.5e4, 3.0e4];
    [0  1.0e4  2.0e4  3.0e4  4.0e4];
    [0 1.0e4 2.0e4 3.0e4 4.0e4]};

T = OC_time_ion(num);

FV_ion = dfana.calcF(JV_solutions_ion{num}, 'whole');
FV_el = dfana.calcF(JV_solutions_el{num}, 'whole');
x_values = JV_solutions_ion{num}.x * 1e7;

hold on

% Do background colours
% Arguments of patch give the coordinates of the corners of the polygon to
% be shaded
patch('XData',[0 x_values(500) x_values(500) 0], 'YData',[1e10 1e10 -1e10 -1e10], 'FaceColor', colors_JV{1},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
patch('XData',[x_values(500) x_values(600) x_values(600) x_values(500)], 'YData',[1e10 1e10 -1e10 -1e10], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1100) x_values(1200) x_values(1200) x_values(1100)], 'YData',[1e10 1e10 -1e10 -1e10], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1200) x_values(end) x_values(end) x_values(1200)], 'YData',[1e10 1e10 -1e10 -1e10],  'FaceColor', colors_JV{3},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')

plot(x_values, -FV_ion(T,:), 'Color', [0.4940 0.1840 0.5560])
plot(x_values, -FV_el(T,:), 'Color', [0.4940 0.1840 0.5560], 'LineStyle', '--')
plot(x_values, zeros(length(x_values)), 'Color', 'black', 'LineWidth', 2)

hold off

box on
xlim([0, max(x_values)])
ylim([lower_lim(num), upper_lim(num)])
xticks([0 100 200 300 400])
%yticks(ytick_cell{num})
ylabel('|Electric Field Strength| (Vcm^{-1})', 'FontSize', fontsize)
xlabel('Perovskite depth (nm)', 'FontSize', fontsize)
legend('', '', '', '', 'F_{ion}', 'F_{no ions}', '', 'NumColumns', 1, 'Position', [0.68 0.8 0.15 0.02], 'FontSize', fontsize)

fig1c = gcf;


%% Plot carrier populations at SC and OC, with and without mobile ions
% and the ion population at OC 

num = 2;
x_values = JV_solutions_ion{num}.x * 1e7;
if num ~= 1
    pos_el = [0.68 0.8 0.1 0.05];
    pos_ion = [0.66 0.8 0.1 0.05];
else
    pos_el = [0.7 0.8 0.1 0.05];
    pos_ion = [0.68 0.8 0.1 0.05];
end

%SC
T = 31;

figure('Name', 'Carrier Densities SC', 'Position', [100 100 1250 2000])
set(gca, 'DefaultLineLineWidth', lw, 'FontSize', ticksize)

hold on

% Do background colours
% Arguments of patch give the coordinates of the corners of the polygon to
% be shaded
patch('XData',[0 x_values(500) x_values(500) 0], 'YData',[30 30 2 2], 'FaceColor', colors_JV{1},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
patch('XData',[x_values(500) x_values(600) x_values(600) x_values(500)], 'YData',[30 30 2 2], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1100) x_values(1200) x_values(1200) x_values(1100)], 'YData',[30 30 2 2], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1200) x_values(end) x_values(end) x_values(1200)], 'YData',[30 30 2 2],  'FaceColor', colors_JV{3},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')

semilogy(x_values, log10(JV_solutions_ion{num}.u(T,:,2)), 'b-',...
    x_values, log10(JV_solutions_ion{num}.u(T,:,3)), 'r-',...    
    x_values, log10(JV_solutions_el{num}.u(T,:,2)), 'b--',...
    x_values, log10(JV_solutions_el{num}.u(T,:,3)), 'r--')

ylim([9,20])
if num == 1
    xlim([0, x_values(1540)])
else
    xlim([0, max(x_values)])
end
ylabel('log_{10}(Carrier Density/cm^{-3})', 'Fontsize', fontsize)
xlabel('Device depth (nm)', 'Fontsize', fontsize)
yticks([10 12 14 16 18 20])
xticks([0 100 200 300 400])
legend('','','','','electrons','holes', 'NumColumns', 1, 'Position', pos_el, 'Fontsize', fontsize)
title(legend, 'Short Circuit', 'Fontsize', fontsize)

hold off 

fig2 = gcf;

%OC (el)
T = OC_time_el(num);

figure('Name', 'Carrier Densities OC', 'Position', [100 100 1250 2000])
set(gca, 'DefaultLineLineWidth', lw, 'FontSize', ticksize)

hold on

% Do background colours
% Arguments of patch give the coordinates of the corners of the polygon to
% be shaded
patch('XData',[0 x_values(500) x_values(500) 0], 'YData',[30 30 2 2], 'FaceColor', colors_JV{1},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
patch('XData',[x_values(500) x_values(600) x_values(600) x_values(500)], 'YData',[30 30 2 2], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1100) x_values(1200) x_values(1200) x_values(1100)], 'YData',[30 30 2 2], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1200) x_values(end) x_values(end) x_values(1200)], 'YData',[30 30 2 2],  'FaceColor', colors_JV{3},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')

semilogy(x_values, log10(JV_solutions_ion{num}.u(T,:,2)), 'b-',...
         x_values, log10(JV_solutions_ion{num}.u(T,:,3)), 'r-')
semilogy(x_values, log10(JV_solutions_el{num}.u(T,:,2)), 'b--',...
         x_values, log10(JV_solutions_el{num}.u(T,:,3)), 'r--')

hold off 

ylim([9,20])
if num == 1
    xlim([0, x_values(1540)])
else
    xlim([0, max(x_values)])
end
ylabel('log_{10}(Carrier Density/cm^{-3})', 'Fontsize', fontsize)
xlabel('Device depth (nm)', 'Fontsize', fontsize)
yticks([10 12 14 16 18 20])
xticks([0 100 200 300 400])
legend('','','','','electrons','holes', 'NumColumns', 1, 'Position', pos_el, 'Fontsize', fontsize)
if num == 1
    leg_title = [num2str(Voc_el(num), 2) ' V'];
else
    leg_title = [num2str(Voc_el(num), 3) ' V'];
end
title(legend, leg_title, 'Fontsize', fontsize)

fig3a = gcf;

%OC (ion)
T = OC_time_ion(num);

figure('Name', 'Carrier Densities OC', 'Position', [100 100 1250 2000])
set(gca, 'DefaultLineLineWidth', lw, 'FontSize', ticksize)

hold on

% Do background colours
% Arguments of patch give the coordinates of the corners of the polygon to
% be shaded
patch('XData',[0 x_values(500) x_values(500) 0], 'YData',[30 30 2 2], 'FaceColor', colors_JV{1},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
patch('XData',[x_values(500) x_values(600) x_values(600) x_values(500)], 'YData',[30 30 2 2], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1100) x_values(1200) x_values(1200) x_values(1100)], 'YData',[30 30 2 2], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1200) x_values(end) x_values(end) x_values(1200)], 'YData',[30 30 2 2],  'FaceColor', colors_JV{3},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')

semilogy(x_values, log10(JV_solutions_ion{num}.u(T,:,2)), 'b-',...
         x_values, log10(JV_solutions_ion{num}.u(T,:,3)), 'r-')
semilogy(x_values, log10(JV_solutions_el{num}.u(T,:,2)), 'b--',...
         x_values, log10(JV_solutions_el{num}.u(T,:,3)), 'r--')

hold off 

ylim([9,20])
if num == 1
    xlim([0, x_values(1540)])
else
    xlim([0, max(x_values)])
end
ylabel('log_{10}(Carrier Density/cm^{-3})', 'Fontsize', fontsize)
xlabel('Device depth (nm)', 'Fontsize', fontsize)
yticks([10 12 14 16 18 20])
xticks([0 100 200 300 400])
legend('','','','','electrons','holes', 'NumColumns', 1, 'Position', pos_el, 'Fontsize', fontsize)
leg_title = [num2str(Voc_ion(num), 3) ' V'];
title(legend, leg_title, 'Fontsize', fontsize)

fig3b = gcf;


%Ions OC
figure('Name', 'Cation Densities OC', 'Position', [100 100 1250 2000])
set(gca, 'DefaultLineLineWidth', lw, 'FontSize', ticksize)

hold on

% Do background colours
% Arguments of patch give the coordinates of the corners of the polygon to
% be shaded
patch('XData',[0 x_values(500) x_values(500) 0], 'YData',[1e19 1e19 1e9 1e9], 'FaceColor', colors_JV{1},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
patch('XData',[x_values(500) x_values(600) x_values(600) x_values(500)], 'YData',[1e19 1e19 1e9 1e9], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1100) x_values(1200) x_values(1200) x_values(1100)], 'YData',[1e19 1e19 1e9 1e9], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1200) x_values(end) x_values(end) x_values(1200)], 'YData',[1e19 1e19 1e9 1e9],  'FaceColor', colors_JV{3},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')

T = OC_time_ion(num);

semilogy(x_values, (JV_solutions_ion{num}.u(T,:,4)), 'color', colors_JV{2}, 'LineStyle', '-')

hold off 

ylim([1e8, 1e18])
if num == 1
    xlim([0, x_values(1540)])
else
    xlim([0, max(x_values)])
end
ylabel('Cation Density (cm^{-3})', 'Fontsize', fontsize)
xlabel('Device depth (nm)', 'Fontsize', fontsize)
xticks([0 100 200 300 400])
yticks([2e17 4e17 6e17 8e17 10e17])
legend('','','','','Open Circuit', 'NumColumns', 1, 'Position', pos_ion, 'Fontsize', fontsize)
title(legend, 'Ion Distribution', 'Fontsize', fontsize)

fig4 = gcf;

%% Save Images
%Set details of what you're saving
fig_num = 1.5;
ETM_distplot = 2;
ETM_fieldplot = 2;

if save == 1 && fig_num == 1
    filename = ['C:\Users\ljh3218\OneDrive - Imperial College London\PhD\IonEfficiency\ESAFigures_1Sun\FieldScreening_ETM' num2str(ETM_fieldplot) '.png'];
    exportgraphics(fig1a, filename, 'Resolution', 300)
elseif save == 1 && fig_num == 1.5
    filename = ['C:\Users\ljh3218\OneDrive - Imperial College London\PhD\IonEfficiency\ESAFigures_1Sun\CationDistributionsSC_ETM' num2str(ETM_distplot) '.png'];
    exportgraphics(fig1b, filename, 'Resolution', 300)
elseif save == 1 && fig_num == 2
    filename = ['C:\Users\ljh3218\OneDrive - Imperial College London\PhD\IonEfficiency\ESAFigures_1Sun\CarrierDistributionsSC_ETM' num2str(ETM_distplot) '.png'];
    exportgraphics(fig2, filename, 'Resolution', 300)
elseif save == 1 && fig_num == 3
    filename = ['C:\Users\ljh3218\OneDrive - Imperial College London\PhD\IonEfficiency\ESAFigures_1Sun\CarrierDistributionsOCel_ETM' num2str(ETM_distplot) '.png'];
    exportgraphics(fig3a, filename, 'Resolution', 300)
elseif save == 1 && fig_num == 3.5
    filename = ['C:\Users\ljh3218\OneDrive - Imperial College London\PhD\IonEfficiency\ESAFigures_1Sun\CarrierDistributionsOCion_ETM' num2str(ETM_distplot) '.png'];
    exportgraphics(fig3b, filename, 'Resolution', 300)
elseif save == 1 && fig_num == 4
    filename = ['C:\Users\ljh3218\OneDrive - Imperial College London\PhD\IonEfficiency\ESAFigures_1Sun\CationDistributionsOC_ETM' num2str(ETM_distplot) '.png'];
    exportgraphics(fig4, filename, 'Resolution', 300)
end