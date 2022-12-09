%Code to model the results of Weidong's experiment varying the choice of
%ETL material used in his FaCs devices.

%% Read in data files 
par1 = pc('Input_files/HTL_MAPI_NoOffset.csv');
par2 = pc('Input_files/HTL_MAPI_NegOffset.csv');
par3 = pc('Input_files/HTL_MAPI_PosOffset.csv');

devices = {par1, par2, par3};
num_devices = length(devices);

%% Find eqm solutions 
eqm_solutions_dark = cell(1,num_devices);
for i = 1:num_devices
    eqm_solutions_dark{i} = equilibrate(devices{i});
end

%% Perform CV scans
CV_solutions_ion = cell(1,num_devices);
for j = 1:num_devices
    suns = 1;
    CV_solutions_ion{j} = doCV(eqm_solutions_dark{j}.ion, suns, -0.3, 1.3, -0.3, 10e-3, 1, 321);
end

%% Plot JVs
figure('Name', 'JVPlot', 'Position', [100 100 1250 2000])
%Colour order is red, green, blue
colors_JV = {[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],[0 0.4470 0.7410]};

hold on
xline(0, 'color', 'black', 'LineWidth', 1)
yline(0, 'color', 'black', 'LineWidth', 1)

for m=1:3
    order = [2 1 3];
    v = dfana.calcVapp(CV_solutions_ion{order(m)});
    j = dfana.calcJ(CV_solutions_ion{order(m)}).tot(:,1);
    
    plot(v(:), j(:)*1000, 'color', colors_JV{m}, 'LineWidth', 3) 
end   

hold off

set(gca, 'FontSize', 25)
xlim([-0.25, 1.25])
ylim([-23, 5])
xlabel('Voltage (V)', 'FontSize', 30)
ylabel('Current Density (mAcm^{-2})', 'FontSize', 30)
title(legend, 'E_{LUMO} - E_{CB}')
legend({'', '', ' -0.2 eV', '  0.0 eV', '+0.2 eV'}, 'Location', 'northwest', 'FontSize', 30)
fig1 = gcf;
%% Break down contributions to the current
%Columns in J_values are J_gen, J_rad, J_srh, J_vsr (left), J_vsr (right) and J_ext

num_values = length(CV_solutions_ion{1}.t);
J_values = zeros(num_values,7,num_devices);
J_values_el = zeros(num_values,7,num_devices);
Voc_values_ion = zeros(num_devices);
e = -CV_solutions_ion{1}.par.e;

for k=1:num_devices
    CVsol = CV_solutions_ion{k};
    loss_currents = dfana.calcr(CVsol,'sub');
    x = CVsol.par.x_sub;
    num_points = length(x);
    gxt = dfana.calcg(CVsol);
    J = dfana.calcJ(CVsol);
    j_surf_rec = dfana.calcj_surf_rec(CVsol);

    J_values(:,1,k) = e*trapz(x, gxt(1,:))';
    J_values(:,2,k) = e*trapz(x, loss_currents.btb, 2)';
    J_values(:,3,k) = e*trapz(x, loss_currents.srh, 2)';
    J_values(:,4,k) = e*trapz(x(1:ceil(num_points/2)), loss_currents.vsr(:,1:ceil(num_points/2)), 2)';
    J_values(:,5,k) = e*trapz(x(ceil(num_points/2)+1:end), loss_currents.vsr(:,ceil(num_points/2)+1:end), 2)';
    J_values(:,6,k) = e*(j_surf_rec.tot);
    J_values(:,7,k) = J.tot(:,1);

    Voc_values_ion(k) = CVstats(CVsol).Voc_f;

end    

%% Plot contributons to the current
%J_rad not corrected for EL 
figure('Name', 'Current Contributions')
num = 3;
line_colour = {[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560]...
                [0 0.4470 0.7410], [0.3010 0.7450 0.9330], 'black', [0.4660 0.6740 0.1880]};
V = dfana.calcVapp(CV_solutions_ion{1});

box on
hold on
for n = 1:7
    plot(V(:), J_values(:,n,num), 'color', line_colour{n}) 
end
plot(V(1:num_values), zeros(1,num_values), 'black', 'LineWidth', 1)
xline(0, 'color', 'black', 'LineWidth', 1)
yline(0, 'color', 'black', 'LineWidth', 1)
hold off
xlim([-0.25, 1.25])
xlabel('Voltage (V)')
ylim([-0.023, 0.005])
ylabel('Current Density (Acm^{-2})')
legend({'J_{gen}', 'J_{rad}', 'J_{SRH}', 'J_{interface (left)}', 'J_{interface (right)}', '','J_{ext}'}, 'Location', 'bestoutside')

fig2 = gcf;

%% Plot Jnonrad/Jrad
%J_rad not corrected for EL - see EL_Measurements
figure('Name', 'Jnonrad/Jrad', 'Position', [100 100 1250 2000])

box on
for n = 1:num_devices
    order = [2 1 3];
    semilogy(V(1:ceil(num_values/2)), sum(J_values(1:ceil(num_values/2),3:6,order(n)),2)./J_values(1:ceil(num_values/2),2,order(n)), 'color', colors_JV{n}, 'LineWidth', 3) 
    hold on
end
xline(0, 'color', 'black', 'LineWidth', 1)
yline(0, 'color', 'black', 'LineWidth', 1)
hold off
set(gca, 'FontSize', 25)
xlim([-0.25, 1.25])
xlabel('Voltage (V)', 'FontSize', 30)
ylim([5, 1200])
ylabel('J_{non rad}/J_{rad}', 'FontSize', 30)
legend({'  - 0.2 eV', '    0.0 eV', ' + 0.2 eV'}, 'Location', 'southwest', 'FontSize', 30)

fig3 = gcf;
%% Save Plots at 300 dpi
save = 1;
fig_num = 3;

if save == 1 && fig_num == 1
    exportgraphics(fig1, ...
    'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Paper1p5\Simulations\Ideal_HTL\JV.png', ...
    'Resolution', 300)
elseif save == 1 && fig_num == 2
    filename = ["C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Paper1p5\Simulations\Ideal_HTL\Current_Contribution_ETL" + num2str(num) + ".png"];
    exportgraphics(fig2, filename, 'Resolution', 300)
elseif save == 1 && fig_num == 3
    exportgraphics(fig3, ...
    'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Paper1p5\Simulations\Ideal_HTL\Jnonrad_Jrad.png', ...
    'Resolution', 300)
end 










