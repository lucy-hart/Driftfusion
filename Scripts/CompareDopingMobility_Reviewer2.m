%Code to model the results of Weidong's experiment varying the choice of
%ETL material used in his FaCs devices.

%% Read in data files 
par_pcbm = pc('Input_files/PTAA_MAPI_PCBM_v5.csv');
par_pcbm_LowestMobility = pc('Input_files/PTAA_MAPI_PCBM_LowestMob.csv');
par_pcbm_LowestMobility_Doped = pc('Input_files/PTAA_MAPI_PCBM_DopedETM.csv');

devices = {par_pcbm, par_pcbm_LowestMobility, par_pcbm_LowestMobility_Doped};
num_devices = length(devices);

%Change light source to laser 
for i = 1:num_devices
    par = devices{i};
    par.light_source1 = 'laser';
    par.laser_lambda1 = 532;
    par.pulsepow = 62;
    par.RelTol_vsr = 0.1;
    devices{i} = refresh_device(par);
end

%% Find eqm solutions 
eqm_solutions_dark = cell(1,num_devices);
for i = 1:num_devices
    eqm_solutions_dark{i} = equilibrate(devices{i});
end

%% Perform CV scans
CV_solutions_el = cell(1,num_devices);
CV_solutions_ion = cell(1,num_devices);
for j = 1:num_devices
    suns = 1;
    CV_solutions_ion{j} = doCV(eqm_solutions_dark{j}.ion, suns, -0.3, 1.3, -0.3, 10e-3, 1, 321);
end

%% Plot JVs
figure('Name', 'JVPlot', 'Position', [100 100 1250 2000])
%Colour order is green, blue, red
colors_JV = {[0.4660 0.6740 0.1880],[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};

for m=1:3
    v = dfana.calcVapp(CV_solutions_ion{m});
    j = -dfana.calcJ(CV_solutions_ion{m}).tot(:,1);
    plot(v(:), j(:)*1000, 'color', colors_JV{m}, 'LineWidth', 3) 
    hold on
end
plot(v(:), zeros(1,length(v)), 'black', 'LineWidth', 1)
hold off

set(gca, 'FontSize', 35)
xlim([0, 1.3])
ylim([0,27])
xlabel('Voltage (V)', 'FontSize', 40)
ylabel('Current Density (mAcm^{-2})', 'FontSize', 40)
legend({'  ETM 1', '  ETM 2', '  ETM 3'}, 'Location', 'southwest', 'FontSize', 40)
ax1 = gca;

%% Break down contributions to the current
%Columns in J_values are J_gen, J_rad, J_srh, J_vsr (left), J_vsr (right) and J_ext

num_values = length(CV_solutions_ion{1}.t);
J_values = zeros(num_values,7,num_devices);
Voc_values_ion = zeros(num_devices,2);
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
    J_values(:,4,k) = e*trapz(x(1:num_points), loss_currents.vsr(:,1:num_points), 2)';
    J_values(:,5,k) = e*trapz(x(num_points+1:end), loss_currents.vsr(:,num_points+1:end), 2)';
    J_values(:,6,k) = e*(j_surf_rec.tot);
    J_values(:,7,k) = J.tot(:,1);

    Voc_values_ion(k,1) = CVstats(CVsol).Voc_r;
    Voc_values_ion(k,2) = interp1(v(ceil(num_values/2):end), J_values(ceil(num_values/2):end,2,k), Voc_values_ion(k,1));

end      

%% Plot PLQY results
%One plot for effect of changing mobility 
%Another plot for effect of changing energetics
figure('Name', 'PLQYPlot Mobility', 'Position', [100 100 1250 2000])
mob_order = [2,3,1];
V = dfana.calcVapp(CV_solutions_ion{1});
for i = 1:3
    semilogy(V(161:end), 100*(J_values(161:end,2,i))./J_values(161:end,1,i), 'color', colors_JV{i}, 'LineWidth', 3) 
    hold on
end
hold off

set(gca, 'FontSize', 25)
xlim([0, 1.2])
ylim([1e-5, 0.3])
legend({'  10^{-3}', '  10^{-5}', '  10^{-5}, doped'}, 'Location', 'southeast', 'FontSize', 30)
xlabel('Voltage (V)', 'FontSize', 30)
ylabel('PLQY (%)', 'FontSize', 30)
title(legend, 'ETM Mobility (cm^{2}/V/s)', 'Fontsize', 30)
ax2 = gca;

%% Save Plots at 300 dpi
save = 0;
fig_num = 1;

if save == 1 && fig_num == 1
    exportgraphics(ax1, ...
    'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Weidong_ETL\Paper\Reviewer2_HigherVoc.png', ...
    'Resolution', 300)
end
