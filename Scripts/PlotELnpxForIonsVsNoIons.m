% Plot energy levels and carrier popultions at short circuit 
% with and without ion motion.
% SC conditions for forward scan are at 31st timepoint, and 291st for reverse 

num = 3;

[Ecb_ion, Evb_ion, Efn_ion, Efp_ion] = dfana.calcEnergies(CV_solutions_ion{num});
[Ecb_el, Evb_el, Efn_el, Efp_el] = dfana.calcEnergies(CV_solutions_el{num});
x_values = CV_solutions_ion{num}.x * 1e7;

figure(5)
subplot(2,2,1)
plot(x_values, Ecb_ion(31,:), 'b',...
    x_values, Evb_ion(31,:), 'r',...
    x_values, Efn_ion(31,:), 'b--', ...
    x_values, Efp_ion(31,:), 'r--')
xlim([0, max(x_values)])
ylim([-7.5, -2.5])
ylabel('Energy, eV')
xlabel('Device depth (nm)')
title('Energy Levels with Ions')
subplot(2,2,2)
plot(x_values, Ecb_el(31,:), 'b',...
    x_values, Evb_el(31,:), 'r',...
    x_values, Efn_el(31,:), 'b--', ...
    x_values, Efp_el(31,:), 'r--')
xlim([0, max(x_values)])
ylim([-7.5, -2.5])
ylabel('Energy, eV')
xlabel('Device depth (nm)')
legend('E_{C}', 'E_{V}', 'E_{f,n}', 'E_{f,p}', 'NumColumns', 2, 'Location', 'southwest')
title('Energy levels without Ions')
subplot(2,2,[3,4])
semilogy(x_values, CV_solutions_ion{num}.u(31,:,2), 'b',...
    x_values, CV_solutions_ion{num}.u(31,:,3), 'r',...
    x_values, CV_solutions_el{num}.u(31,:,2), 'b--',...
    x_values, CV_solutions_el{num}.u(31,:,3), 'r--')
xlim([0, max(x_values)])
ylim([1e8,1e17])
ylabel('Carrier concentrations (cm^{-3})')
xlabel('Device depth (nm)')
legend('n, ions', 'p, ions', 'n, no ions', 'p, no ions', 'Location', 'bestoutside')


