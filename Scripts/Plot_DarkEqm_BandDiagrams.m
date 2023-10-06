%% File to plot dark eqm energy levels for Davide dark currents
ion = 1;

if ion == 1
    devices = {eqm_PM6.ion, eqm_PBDBT.ion};
elseif ion == 0
    devices = {eqm_PM6.el, eqm_PBDBT.el};
end
energies = cell(4, length(devices));

for i = 1:length(devices)
    [energies{1,i},energies{2,i},energies{3,i},energies{4,i}] = dfana.calcEnergies(devices{i});
end 

figure('Name', 'Dark Eqm Energy Levels')
colours = {[0 0 1], [1 0 0], [0 0 0], [0 0 0]};

hold on
for i = 1:2
    plot(1e7*devices{1}.x, energies{i,1}(end,:), 'Color', colours{i})
    plot(1e7*devices{1}.x, energies{i,2}(end,:), 'Color', colours{i}, 'LineStyle', '--')
end
yline(-5.2, 'Color', 'Black', 'LineWidth', 2)
hold off

xlabel('Distance (nm)')
xlim([0, 1e7*max(devices{1}.x)])
ylabel('Energy (eV)')