%% Ion Density Sweep for Weidong Paper
par = pc('Input_files/PTAA_MAPI_PCBM_v4.csv');
eqm_QJV = equilibrate(par);

ion_densities = [1e16, 5e16, 1e17, 5e17, 1e18];
solutions = cell(1,5);

for m = 1:length(ion_densities)

    par.Ncat(:) = ion_densities(m);
    par.Nani(:) = ion_densities(m);
    par.RelTol_vsr = 0.1;
    par = refresh_device(par);
    eqm_QJV = equilibrate(par);
    solutions{m} = doCV(eqm_QJV.ion, 1.15, -0.3, 1.2, -0.3, 10e-3, 1, 321);

end

%%
figure('Name', 'IonDensityPLQY', 'Position', [100 100 1250 2000])
%Colour order is red, blue, green, yellow, purple
colors_JV = {[0.8500 0.3250 0.0980],[0 0.4470 0.7410],[0.4660 0.6740 0.1880],[0.4940 0.1840 0.5560],[0.9290 0.6940 0.1250]};

box on 
for m = 1:length(ion_densities)

    v = dfana.calcVapp(solutions{m});
    x = solutions{m}.par.x_sub;
    gxt = dfana.calcg(solutions{m});
    J_gen = trapz(x, gxt(1,:))';
    loss_currents = dfana.calcr(solutions{m},'sub');
    J_rad = trapz(x, loss_currents.btb, 2)';

    semilogy(v(:), 100*J_rad(:)./J_gen(:), 'color', colors_JV{m}, 'LineWidth', 3)
    hold on

end
hold off

set(gca, 'FontSize', 25)
xlim([0, 1.2])
ylim([5e-5, 1])
xlabel('Voltage (V)', 'FontSize', 30)
ylabel('PLQY (%)', 'FontSize', 30)
title(legend, 'Ion Density (cm^{-3})', 'Fontsize', 30)
legend({'  1 x 10^{16}', '  5 x 10^{16}', '  1 x 10^{17}', '  5 x 10^{17}', '  1 x 10^{18}'}, 'Location', 'northwest', 'FontSize', 30)
ax1 = gca;

%% Save Plot
save = 1;

if save == 1 
    exportgraphics(ax1, ...
    'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Weidong_ETL\Paper\v2\SweepIonDensity.png', ...
    'Resolution', 300)
end 