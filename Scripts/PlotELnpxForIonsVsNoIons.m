%Need to run this after CalculateQFLS so that OC_time has been calculated 

%Plot parameters
lw = 3;
fontsize = 30;
ticksize = 20;
save = 1;
%% Plot energy levels and carrier popultions at open circuit for PS1 and PS2
% to illustrate band bending at perovskite/ETM interface cause by deep LUMO
% of ETM in PS1. 

%set x limits
numstart = 1050;
numstop_PS1 = 1325;
numstop_PS2 = 1300;

%PS1
num = 1;
T = OC_time_ion(num);
[Ecb, Evb, Efn, Efp] = dfana.calcEnergies(CV_solutions_ion{num});
x_values = CV_solutions_ion{num}.x * 1e7;

%Make band bending more obvious
%Shift value up by 0.2 eV in the bulk and modulate by a sigmoid
%This is bad science
Efn_dramatic_bending = Efn(T,numstart:1100) + 0.2./(1+exp((x_values(numstart:1100)-411.9)/0.02));

figure('Name', 'Band Bending Kloc6 Inset', 'Position', [100 100 1250 2000])
set(gca, 'DefaultLineLineWidth', 5)
num_crop = 25;

hold on

% Do background colours
% Arguments of patch give the coordinates of the corners of the polygon to
% be shaded
patch('XData',[x_values(1200)-2 x_values(end)-2 x_values(end)-2 x_values(1200)-2], ...
    'YData',[-1 -1 -10 -10],  'FaceColor', colors_JV{3},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')

plot(x_values(numstart+num_crop:1100), Ecb(T,numstart+num_crop:1100), 'b',...
    x_values(1200:numstop_PS1-num_crop)-2, Ecb(T,1200:numstop_PS1-num_crop), 'b',...
    x_values(numstart+num_crop:1100), Efn(T,numstart+num_crop:1100), 'b--',...
    x_values(1200:numstop_PS1-num_crop)-2, Efn(T,1200:numstop_PS1-num_crop), 'b--')
plot(x_values(numstart+num_crop:1100-1), Efn(T,numstart).*ones(1,length(x_values(numstart+num_crop:1100))-1), ...
    'color', [0.9290 0.6940 0.1250], 'LineStyle', '--')

hold off

xlim([x_values(numstart+num_crop), x_values(numstop_PS1-num_crop)-2])
ylim([-4.2, -3.6])
set(gca, 'YTick',[], 'XTick', [])
box on
leg = legend('','', '', '', '', 'E_{f,n} (bulk)', ...
    'NumColumns', 1, 'Location', 'west', 'Fontsize', 50);
leg.ItemTokenSize = [80, 18];
fig1a_inset = gcf;

figure('Name', 'Band Bending Kloc6', 'Position', [100 100 1250 2000])
set(gca, 'DefaultLineLineWidth', lw)

hold on

% Do background colours
% Arguments of patch give the coordinates of the corners of the polygon to
% be shaded
patch('XData',[x_values(1200)-2 x_values(end)-2 x_values(end)-2 x_values(1200)-2], ...
    'YData',[-1 -1 -10 -10],  'FaceColor', colors_JV{3},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')

plot(x_values(numstart:1100), Ecb(T,numstart:1100), 'b',...
    x_values(1200:numstop_PS1)-2, Ecb(T,1200:numstop_PS1), 'b',...
    x_values(numstart:1100), Evb(T,numstart:1100), 'r',...
    x_values(1200:numstop_PS1)-2, Evb(T,1200:numstop_PS1), 'r',...
    x_values(numstart:1100), Efn(T,numstart:1100), 'b--',...
    x_values(1200:numstop_PS1)-2, Efn(T,1200:numstop_PS1), 'b--',...
    x_values(numstart:1100), Efp(T,numstart:1100), 'r--',...
    x_values(1200:numstop_PS1)-2, Efp(T,1200:numstop_PS1), 'r--')

hold off

xlim([x_values(numstart), x_values(numstop_PS1)-2])
ylim([-6.5, -3.5])
set(gca, 'YTick',[], 'XTick', [])
box on
ylabel('Energy', 'Fontsize', fontsize)
xlabel('Device depth', 'Fontsize', fontsize)
legend('','E_{C}', '', 'E_{V}', '', 'E_{f,n}', '', 'E_{f,p}', '', ...
    'NumColumns', 2, 'Location', 'southwest', 'Fontsize', fontsize)

fig1a = gcf;

%PS2 
num = 2;
T = OC_time_ion(num);
[Ecb_ion, Evb_ion, Efn_ion, Efp_ion] = dfana.calcEnergies(CV_solutions_ion{num});
FV = dfana.calcF(CV_solutions_ion{num}, 'whole');
x_values = CV_solutions_ion{num}.x * 1e7;

figure('Name', 'Band Bending PCBM', 'Position', [100 100 1250 2000])
set(gca, 'DefaultLineLineWidth', lw)

hold on

% Do background colours
% Arguments of patch give the coordinates of the corners of the polygon to
% be shaded
patch('XData',[x_values(1200)-2 x_values(end)-2 x_values(end)-2 x_values(1200)-2], ...
    'YData',[-1 -1 -10 -10],  'FaceColor', colors_JV{3},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')

plot(x_values(numstart:1100), Ecb_ion(T,numstart:1100), 'b',...
    x_values(1200:numstop_PS2)-2, Ecb_ion(T,1200:numstop_PS2), 'b',...
    x_values(numstart:1100), Evb_ion(T,numstart:1100), 'r',...
    x_values(1200:numstop_PS2)-2, Evb_ion(T,1200:numstop_PS2), 'r',...
    x_values(numstart:1100), Efn_ion(T,numstart:1100), 'b--',...
    x_values(1200:numstop_PS2)-2, Efn_ion(T,1200:numstop_PS2), 'b--',...
    x_values(numstart:1100), Efp_ion(T,numstart:1100), 'r--',...
    x_values(1200:numstop_PS2)-2, Efp_ion(T,1200:numstop_PS2), 'r--')

hold off 

xlim([x_values(numstart), x_values(numstop_PS2)-2])
ylim([-6.5, -3.5])
set(gca, 'YTick',[], 'XTick', [])
box on 
ylabel('Energy', 'Fontsize', fontsize)
xlabel('Device depth', 'Fontsize', fontsize)
legend('','E_{C}', '', 'E_{V}', '', 'E_{f,n}', '', 'E_{f,p}', '', ...
    'NumColumns', 2, 'Location', 'southwest', 'Fontsize', fontsize)

fig1b = gcf;
%% Plot carrier populations at SC (top) and OC (middle), with and without mobile ions
% and the ion population at OC vs SC (bottom)

num = 2;
x_values = CV_solutions_ion{num}.x * 1e7;
if num == 2
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

semilogy(x_values, log10(CV_solutions_ion{num}.u(T,:,2)), 'b-',...
    x_values, log10(CV_solutions_ion{num}.u(T,:,3)), 'r-',...    
    x_values, log10(CV_solutions_el{num}.u(T,:,2)), 'b--',...
    x_values, log10(CV_solutions_el{num}.u(T,:,3)), 'r--')

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

%OC
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

semilogy(x_values, log10(CV_solutions_ion{num}.u(T,:,2)), 'b-',...
         x_values, log10(CV_solutions_ion{num}.u(T,:,3)), 'r-')

T = OC_time_el(num);

semilogy(x_values, log10(CV_solutions_el{num}.u(T,:,2)), 'b--',...
         x_values, log10(CV_solutions_el{num}.u(T,:,3)), 'r--')

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
title(legend, 'Open Circuit', 'Fontsize', fontsize)

fig3 = gcf;

%Ions
upper_lim_ion = [1e18, 1e18, 1.4e18];
T = 31;

figure('Name', 'Cation Densities', 'Position', [100 100 1250 2000])
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

plot(x_values, (CV_solutions_ion{num}.u(T,:,4)), 'color', colors_JV{2})

T = OC_time_ion(num);

semilogy(x_values, (CV_solutions_ion{num}.u(T,:,4)), 'color', colors_JV{2}, 'LineStyle', '--')

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
if num == 3
    yticks([2e17 4e17 6e17 8e17 10e17 12e17 14e17])
else
    yticks([2e17 4e17 6e17 8e17 10e17])
end
legend('','','','','Short Circuit','Open Circuit', 'NumColumns', 1, 'Position', pos_ion, 'Fontsize', fontsize)
title(legend, 'Ion Distribution', 'Fontsize', fontsize)

fig4 = gcf;

%% Plot E field at SC with and without effects of mobile ions

figure('Name', 'Field Screening', 'Position', [100 100 1250 2000])
set(gca, 'DefaultLineLineWidth', lw, 'FontSize', ticksize)

num = 2;
upper_lim = [2.0e4, 4.5e4, 6.0e4];
ytick_cell = {[0 0.5e4 1.0e4 1.5e4 2.0e4];
    [0  1.0e4  2.0e4  3.0e4  4.0e4];
    [0 1.0e4 2.0e4 3.0e4 4.0e4 5.0e4 6.0e4]};
T = 31;

FV_ion = dfana.calcF(CV_solutions_ion{num}, 'whole');
FV_el = dfana.calcF(CV_solutions_el{num}, 'whole');
x_values = CV_solutions_ion{num}.x * 1e7;

hold on

plot(x_values-12, -FV_ion(T,:), 'Color', [0.4940 0.1840 0.5560])
plot(x_values-12, -FV_el(T,:), 'Color', [0.4940 0.1840 0.5560], 'LineStyle', '--')
plot(x_values-12, zeros(length(x_values)), 'Color', 'black', 'LineWidth', 2)

hold off

xlim([0, 400])
ylim([-0.5e4, upper_lim(num)])
xticks([0 100 200 300 400])
yticks(ytick_cell{num})
ylabel('|Electric Field Strength| (Vcm^{-1})', 'FontSize', fontsize)
xlabel('Perovskite depth (nm)', 'FontSize', fontsize)
legend('F_{ion}', 'F_{no ions}', '', 'NumColumns', 1, 'Position', [0.72 0.8 0.15 0.02], 'FontSize', fontsize)

fig5 = gcf;
%% Save Images
%Set details of what you're saving
fig_num = 1;
PS_distplot = 1;
PS_fieldplot = 2;

if save == 1 && fig_num == 1
    exportgraphics(fig1a, ...
    'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Weidong_ETL\Paper\BandBending_Kloc6.png', ...
    'Resolution', 300)
    exportgraphics(fig1a_inset, ...
    'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Weidong_ETL\Paper\BandBending_Kloc6_inset.png', ...
    'Resolution', 300)
    exportgraphics(fig1b, ...
    'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Weidong_ETL\Paper\BandBending_PCBM.png', ...
    'Resolution', 300)
elseif save == 1 && fig_num == 2
    filename = ['C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Weidong_ETL\Paper\CarrierDistributionsSC_PS' num2str(PS_distplot) '.png'];
    exportgraphics(fig2, filename, 'Resolution', 300)
elseif save == 1 && fig_num == 3
    filename = ['C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Weidong_ETL\Paper\CarrierDistributions0C_PS' num2str(PS_distplot) '.png'];
    exportgraphics(fig3, filename, 'Resolution', 300)
elseif save == 1 && fig_num == 4
    filename = ['C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Weidong_ETL\Paper\CationDistributions_PS' num2str(PS_distplot) '.png'];
    exportgraphics(fig4, filename, 'Resolution', 300)
elseif save == 1 && fig_num == 5
    filename = ['C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Weidong_ETL\Paper\FieldScreening_PS' num2str(PS_fieldplot) '.png'];
    exportgraphics(fig5, filename, 'Resolution', 300)
end