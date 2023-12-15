%% File to plot dark eqm energy levels for Davide dark currents
eqm = 0;

devices = {Jdark{1}.sol, Jdark{2}.sol, Jdark{3}.sol};

energies_eqm = cell(4, length(devices));
energies_reverse = cell(4, length(devices));

for i = 1:length(devices)
    [energies_reverse{1,i},energies_reverse{2,i},energies_reverse{4,i},energies_reverse{3,i}] = dfana.calcEnergies(devices{i}{1});
    [energies_eqm{1,i},energies_eqm{2,i},energies_eqm{4,i},energies_eqm{3,i}] = dfana.calcEnergies(devices{i}{2});
end 

%%
num = 3;
x_values = 1e7*devices{num}{2}.x;
figure('Name', 'Dark Eqm Energy Levels', 'Position', [50 50 1800 1300])
tiledlayout(1,2)

colours = {[0 0 1], [1 0 0], [0 0 1], [1 0 0]};
colors_JV = {[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],[0 0.4470 0.7410], [0.3010 0.7450 0.9330], [1 1 1]};

nexttile
hold on

patch('XData',[0 x_values(500) x_values(500) 0], 'YData',[0 0 -8 -8], 'FaceColor', colors_JV{1},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
patch('XData',[x_values(500) x_values(600) x_values(600) x_values(500)], 'YData',[0 0 -8 -8], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(600) x_values(1100) x_values(1100) x_values(600)], 'YData',[0 0 -8 -8], 'FaceColor', colors_JV{5},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1100) x_values(1200) x_values(1200) x_values(1100)], 'YData',[0 0 -8 -8], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1200) x_values(2200) x_values(2200) x_values(1200)], 'YData',[0 0 -8 -8],  'FaceColor', colors_JV{3},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
patch('XData',[x_values(2200) x_values(2300) x_values(2300) x_values(2200)], 'YData',[0 0 -8 -8], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(2300) x_values(end) x_values(end) x_values(2300)], 'YData',[0 0 -8 -8], 'FaceColor', colors_JV{4},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')

for i = 1:2
    plot(x_values, energies_eqm{i,num}(end,:), 'Color', colours{i})
end
yline(-5.2, 'Color', 'Black', 'LineWidth', 2, 'LineStyle', '--')
hold off

set(gca, 'Fontsize', 25)
xlabel('Distance (nm)', 'FontSize', 30)
xlim([0, x_values(end)])
ylim([-8, -2])
ylabel('Energy (eV)', 'FontSize', 30)
legend({'', '', '', '', '', '', '','E_{C}', 'E_{V}', 'E_{F}'}, 'FontSize', 25)
title('Dark Equilibrium', 'FontSize', 30)

nexttile
hold on

patch('XData',[0 x_values(500) x_values(500) 0], 'YData',[0 0 -8 -8], 'FaceColor', colors_JV{1},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
patch('XData',[x_values(500) x_values(600) x_values(600) x_values(500)], 'YData',[0 0 -8 -8], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(600) x_values(1100) x_values(1100) x_values(600)], 'YData',[0 0 -8 -8], 'FaceColor', colors_JV{5},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1100) x_values(1200) x_values(1200) x_values(1100)], 'YData',[0 0 -8 -8], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1200) x_values(2200) x_values(2200) x_values(1200)], 'YData',[0 0 -8 -8],  'FaceColor', colors_JV{3},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
patch('XData',[x_values(2200) x_values(2300) x_values(2300) x_values(2200)], 'YData',[0 0 -8 -8], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(2300) x_values(end) x_values(end) x_values(2300)], 'YData',[0 0 -8 -8], 'FaceColor', colors_JV{4},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')

for i = 1:4
    if i <= 2
        plot(x_values, energies_reverse{i,num}(end,:), 'Color', colours{i})
    else
        plot(x_values, energies_reverse{i,num}(end,:), 'Color', colours{i}, 'LineStyle','--')
    end
end
hold off

set(gca, 'Fontsize', 25)
xlabel('Distance (nm)', 'FontSize', 30)
xlim([0, x_values(end)])
ylim([-8, -2])
legend({'', '', '', '', '', '', '', 'E_{C}', 'E_{V}', 'E_{F,n}', 'E_{F,p}'}, 'FontSize', 25)
ylabel('Energy (eV)', 'FontSize', 30)
title('-0.5 V', 'FontSize', 30)

ax1 = gca;
save_fig = 1;
if save_fig == 1
    if num == 1
        exportgraphics(ax1, ...
        'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Davide_OrganicPeroHybrid\FigSxa-DeltaE-0p2eV.png', ...
        'Resolution', 300)
    elseif num == 2
         exportgraphics(ax1, ...
        'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Davide_OrganicPeroHybrid\FigSxa-DeltaE-0p3eV.png', ...
        'Resolution', 300)
    elseif num == 3
         exportgraphics(ax1, ...
        'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Davide_OrganicPeroHybrid\FigSxa-DeltaE-0p4eV.png', ...
        'Resolution', 300)
    end 
end