%% Ion Density Sweep for Weidong Paper
par = pc('Input_files/PTAA_MAPI_PCBM_v5.csv');

ion_densities = [1e16, 5e16, 1e17, 5e17, 1e18];
solutions = cell(1,5);

for m = 1:length(ion_densities)
    
    par.light_source1 = 'laser';
    par.laser_lambda1 = 532;
    par.pulsepow = 62;
    par.RelTol_vsr = 0.1;
    par.Ncat(:) = ion_densities(m);
    par.Nani(:) = ion_densities(m);
    par = refresh_device(par);
    eqm_QJV = equilibrate(par);
    solutions{m} = doCV(eqm_QJV.ion, 1, -0.3, 1.3, -0.3, 10e-3, 1, 321);

end

%%
figure('Name', 'IonDensityPLQY', 'Position', [100 100 1250 2000])
%Colour order is red, blue, green, yellow, purple
colors_JVIonSweep = {[0.8500 0.3250 0.0980],[0 0.4470 0.7410],[0.4660 0.6740 0.1880],[0.4940 0.1840 0.5560],[0.9290 0.6940 0.1250]};

box on 
for m = 1:length(ion_densities)

    v = dfana.calcVapp(solutions{m});
    x = solutions{m}.par.x_sub;
    gxt = dfana.calcg(solutions{m});
    J_gen = trapz(x, gxt(1,:))';
    loss_currents = dfana.calcr(solutions{m},'sub');
    J_rad = trapz(x, loss_currents.btb, 2)';

    semilogy(v(161:end), 100*J_rad(161:end)./J_gen, 'color', colors_JVIonSweep{m}, 'LineWidth', 3)
    hold on

end
hold off

set(gca, 'FontSize', 35)
xlim([0, 1.2])
ylim([5e-5, 1])
xlabel('Voltage (V)', 'FontSize', 40)
ylabel('PLQY (%)', 'FontSize', 40)
title(legend, 'Ion Density (cm^{-3})', 'Fontsize', 40)
legend({'  1 x 10^{16}', '  5 x 10^{16}', '  1 x 10^{17}', '  5 x 10^{17}', '  1 x 10^{18}'}, 'Location', 'northwest', 'FontSize', 40)
ax1 = gca;

%% Calculate and Plot QFLS
num_start = sum(solutions{1}.par.layer_points(1:2))+1;
num_stop = num_start + solutions{1}.par.layer_points(3)-1;
num_values = length(solutions{1}.t);
x = cell(1, length(ion_densities));
for i = 1:length(ion_densities)
    x{i} = solutions{i}.par.x_sub;
end
d = solutions{1}.par.d(3);
QFLS_ion_SweepDensity = zeros(num_values,length(ion_densities));
for y=1:length(ion_densities)
    [Ecb_ion, Evb_ion, Efn_ion, Efp_ion] = dfana.calcEnergies(solutions{y});
    QFLS_ion_SweepDensity(:,y) = trapz(x{y}(num_start:num_stop), Efn_ion(:, num_start:num_stop)-Efp_ion(:,num_start:num_stop),2)/d;
end

% Calculate OC QFLS
Voc_IonSweep = zeros(1,5);
QFLS_OC_IonSweep = zeros(1,length(ion_densities));

for z=1:length(ion_densities)
    Voc_IonSweep(z) = CVstats(solutions{z}).Voc_r;
    QFLS_OC_IonSweep(z) = interp1(v(ceil(num_values/2):end), QFLS_ion(ceil(num_values/2):end,z), Voc_IonSweep(z));
end

figure('Name', 'IonDensityQFLS', 'Position', [100 100 1250 2000])

box on 
for m = 1:length(ion_densities)

    hold on 
    plot(v(161:end), QFLS_ion_SweepDensity(161:end,m), 'color', colors_JVIonSweep{m}, 'LineWidth', 3)
   
end
hold off

set(gca, 'FontSize', 35)
xlim([0, 1.3])
ylim([0.9, 1.3])
xlabel('Voltage (V)', 'FontSize', 40)
ylabel('QFLS (eV)', 'FontSize', 40)
title(legend, 'Ion Density (cm^{-3})', 'Fontsize', 40)
legend({'  1 x 10^{16}', '  5 x 10^{16}', '  1 x 10^{17}', '  5 x 10^{17}', '  1 x 10^{18}'}, 'Location', 'northwest', 'FontSize', 40)
ax2 = gca;
%% Save Plot
save = 1;

if save == 1 
    exportgraphics(ax2, ...
    'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Weidong_ETL\Paper\v2\SweepIonDensity.png', ...
    'Resolution', 300)
end 