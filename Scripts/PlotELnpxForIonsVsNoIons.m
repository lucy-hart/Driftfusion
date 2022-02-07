% Plot energy levels and carrier popultions at short circuit 
% with and without ion motion.
% SC conditions for forward scan are at 31st timepoint, and 291st for reverse 
% OC at 127, 140 and 182 for Kloc-6, PCBM and ICBA respectively

num = 2;
T = 140;

[Ecb_ion, Evb_ion, Efn_ion, Efp_ion] = dfana.calcEnergies(CV_solutions_ion{num});
[Ecb_el, Evb_el, Efn_el, Efp_el] = dfana.calcEnergies(CV_solutions_el{num});
x_values = CV_solutions_ion{num}.x * 1e7;

figure(5)
subplot(2,2,1)
plot(x_values, Ecb_ion(T,:), 'b',...
    x_values, Evb_ion(T,:), 'r',...
    x_values, Efn_ion(T,:), 'b--', ...
    x_values, Efp_ion(T,:), 'r--')
xlim([0, max(x_values)])
ylim([-7.5, -2.0])
ylabel('Energy, eV')
xlabel('Device depth (nm)')
title('Energy Levels with Ions')
subplot(2,2,2)
plot(x_values, Ecb_el(T,:), 'b',...
    x_values, Evb_el(T,:), 'r',...
    x_values, Efn_el(T,:), 'b--', ...
    x_values, Efp_el(T,:), 'r--')
xlim([0, max(x_values)])
ylim([-7.5, -2.0])
ylabel('Energy, eV')
xlabel('Device depth (nm)')
legend('E_{C}', 'E_{V}', 'E_{f,n}', 'E_{f,p}', 'NumColumns', 2, 'Location', 'northeast')
title('Energy levels without Ions')
subplot(2,2,[3,4])
semilogy(x_values, CV_solutions_ion{num}.u(T,:,2), 'b',...
    x_values, CV_solutions_ion{num}.u(T,:,3), 'r',...
    x_values, CV_solutions_el{num}.u(T,:,2), 'b--',...
    x_values, CV_solutions_el{num}.u(T,:,3), 'r--')
xlim([0, max(x_values)])
ylim([1e8,5e17])
ylabel('Carrier concentrations (cm^{-3})')
xlabel('Device depth (nm)')
legend('n, ions', 'p, ions', 'n, no ions', 'p, no ions', 'Location', 'bestoutside')


