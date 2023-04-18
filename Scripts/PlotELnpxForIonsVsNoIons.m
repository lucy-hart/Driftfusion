%Need to run this after CalculateQFLS so that OC_time has been calculated 

%Plot parameters
lw = 3;
fontsize = 40;
ticksize = 35;

%% Plot carrier populations at SC (top) and OC (middle), with and without mobile ions
% and the ion population at OC vs SC (bottom)

num = 1;
x_values = CV_solutions_ion{num}.x * 1e7;
pos_el = [0.62 0.75 0.1 0.05];
pos_ion = [0.66 0.8 0.1 0.05];

%SC
T = 139

figure('Name', 'Carrier Densities SC', 'Position', [100 100 1250 2000])
set(gca, 'DefaultLineLineWidth', lw, 'FontSize', ticksize)

box on
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

semilogy(x_values, log10(CV_solutions_ion{num}.u(T,:,2)), 'b-',...
    x_values, log10(CV_solutions_ion{num}.u(T,:,3)), 'r-',...    
    x_values, log10(CV_solutions_el{num}.u(T,:,2)), 'b--',...
    x_values, log10(CV_solutions_el{num}.u(T,:,3)), 'r--')

ylim([9,20])
xlim([0, max(x_values)])
ylabel('log_{10}(Carrier Density/cm^{-3})', 'Fontsize', fontsize)
xlabel('Device depth (nm)', 'Fontsize', fontsize)
yticks([10 12 14 16 18 20])
xticks([0 100 200 300 400])
legend('','','','','  electrons','  holes', 'NumColumns', 1, 'Position', pos_el, 'Fontsize', fontsize)
%title('ETM Mobility = 10^{-3} cm^{2}/V/s')
title(legend, 'Open Circuit', 'Fontsize', fontsize)

hold off 

fig1 = gcf;

%Ions
upper_lim_ion = 7e17;
T = 31;

figure('Name', 'Cation Densities', 'Position', [100 100 1250 2000])
set(gca, 'DefaultLineLineWidth', lw, 'FontSize', ticksize)

box on
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

plot(x_values, (CV_solutions_ion{num}.u(T,:,4)), 'color', colors_JV{2})

hold off 

ylim([1e8, upper_lim_ion])
xlim([0, max(x_values)])
ylabel('Cation Density (cm^{-3})', 'Fontsize', fontsize)
xlabel('Device depth (nm)', 'Fontsize', fontsize)
xticks([0 100 200 300 400])
yticks([2e17 4e17 6e17])

fig2 = gcf;

%% Plot E field at SC with and without effects of mobile ions

figure('Name', 'Field Screening', 'Position', [100 100 1250 2000])
set(gca, 'DefaultLineLineWidth', lw, 'FontSize', ticksize)

num = 1;
upper_lim = 4.5e4;
ytick_cell = {[0 0.5e4 1.0e4 1.5e4 2.0e4];
    [0  1.0e4  2.0e4  3.0e4  4.0e4];
    [0 1.0e4 2.0e4 3.0e4 4.0e4 5.0e4 6.0e4]};
T = 31;

FV_ion = dfana.calcF(CV_solutions_ion{num}, 'whole');
FV_el = dfana.calcF(CV_solutions_el{num}, 'whole');
x_values = CV_solutions_ion{num}.x * 1e7;

box on
hold on

plot(x_values-10, -FV_ion(T,:), 'Color', [0.4940 0.1840 0.5560])
plot(x_values-10, -FV_el(T,:), 'Color', [0.4940 0.1840 0.5560], 'LineStyle', '--')
plot(x_values-10, zeros(length(x_values)), 'Color', 'black', 'LineWidth', 2)

hold off

xlim([0, 400])
ylim([-0.5e4, upper_lim])
xticks([0 100 200 300 400])
yticks(ytick_cell{2})
ylabel('|Electric Field Strength| (Vcm^{-1})', 'FontSize', fontsize)
xlabel('Perovskite depth (nm)', 'FontSize', fontsize)
legend('  F_{ion}', '  F_{no ions}', '', 'NumColumns', 1, 'Position', [0.70 0.8 0.15 0.02], 'FontSize', fontsize)

fig3 = gcf;

%% Plot recombination currents as function of location in device at Voc
num = 4;
x_values = CV_solutions_ion{num}.x * 1e7;
num_points = length(CV_solutions_ion{num}.t);
pos_ion = [0.66 0.8 0.1 0.05];

figure('Name', 'Recombination vs x', 'Position', [100 100 1250 2000])
set(gca, 'DefaultLineLineWidth', lw, 'FontSize', ticksize)

box on
hold on

% Do background colours
% Arguments of patch give the coordinates of the corners of the polygon to
% be shaded
patch('XData',[0 x_values(500) x_values(500) 0], 'YData',[1e25 1e25 -2e24 -2e24], 'FaceColor', colors_JV{1},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
patch('XData',[x_values(500) x_values(600) x_values(600) x_values(500)], 'YData',[1e25 1e25 -2e24 -2e24], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1100) x_values(1200) x_values(1200) x_values(1100)], 'YData',[1e25 1e25 -2e24 -2e24], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1200) x_values(1700) x_values(1700) x_values(1200)], 'YData',[1e25 1e25 -2e24 -2e24],  'FaceColor', colors_JV{3},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')

rec_currents = dfana.calcr(CV_solutions_ion{num},'sub');
rec_total = dfana.calcj_surf_rec(CV_solutions_ion{num}).tot + rec_currents.btb + rec_currents.vsr + rec_currents.srh;
rec_total_oc = interp1(V(ceil(num_points/2):end), rec_total(ceil(num_points/2):end,:), Voc_values_ion(num,1));

plot(x*1e7, rec_total_oc(1:1700), 'Color', [0.4940 0.1840 0.5560])

xlim([0, 430])
ylim([-2e24, 9e24])
ylabel('Recombination Rate (cm^{-3}s^{-1})', 'Fontsize', fontsize)
xlabel('Device depth (nm)', 'Fontsize', fontsize)
xticks([0 100 200 300 400])
title('E_{LUMO} = -4.15 eV', 'Fontsize', fontsize)

hold off 

fig4 = gcf;

%% Save Images
%Set details of what you're saving
save = 1;
fig_num = 1;

if save == 1 && fig_num == 1
    exportgraphics(fig1, ...
    'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Weidong_ETL\Paper\v2\CarrierDistributionsSC_LowMob.png', ...
    'Resolution', 300)
elseif save == 1 && fig_num == 2
    filename = 'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Weidong_ETL\Paper\v2\CationDistributionSC_LowMob.png';
    exportgraphics(fig2, filename, 'Resolution', 300)
elseif save == 1 && fig_num == 3
    filename = 'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Weidong_ETL\Paper\v2\FieldDistriution_LowestMob.png';
    exportgraphics(fig3, filename, 'Resolution', 300)
end