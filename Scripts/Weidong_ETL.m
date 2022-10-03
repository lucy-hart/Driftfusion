%Code to model the results of Weidong's experiment varying the choice of
%ETL material used in his FaCs devices.

%% Read in data files 
par_pcbm = pc('Input_files/PTAA_MAPI_PCBM_v4.csv');
par_pcbm_LowerMob = pc('Input_files/PTAA_MAPI_PCBM_LowerMob.csv');
par_pcbm_LowestMob = pc('Input_files/PTAA_MAPI_PCBM_LowestMob.csv');
par_pcbm_HigherLUMO = pc('Input_files/PTAA_MAPI_PCBM_HigherLUMO.csv');
par_pcbm_LowerLUMO = pc('Input_files/PTAA_MAPI_PCBM_LowerLUMO.csv');

devices = {par_pcbm, par_pcbm_LowerMob, par_pcbm_LowestMob, par_pcbm_LowerLUMO, par_pcbm_HigherLUMO};
num_devices = length(devices);

%% Find eqm solutions 
eqm_solutions_dark = cell(1,num_devices);
for i = 1:num_devices
    eqm_solutions_dark{i} = equilibrate(devices{i});
end

%% Perform CV scans
CV_solutions_el = cell(1,num_devices);
CV_solutions_ion = cell(1,num_devices);
for j = 1:num_devices
    suns = 1.15;
    CV_solutions_el{j} = doCV(eqm_solutions_dark{j}.el, suns, -0.3, 1.2, -0.3, 10e-3, 1, 301);
    CV_solutions_ion{j} = doCV(eqm_solutions_dark{j}.ion, suns, -0.3, 1.2, -0.3, 10e-3, 1, 301);
end

%% Plot JVs
figure('Name', 'JVPlot', 'Position', [100 100 1250 2000])
%Colour order is red, green, blue
colors_JV = {[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],[0 0.4470 0.7410]};

m=1;
v = dfana.calcVapp(CV_solutions_ion{m});
v_el = dfana.calcVapp(CV_solutions_el{m});
j = -dfana.calcJ(CV_solutions_ion{m}).tot(:,1);
j_el = -dfana.calcJ(CV_solutions_el{m}).tot(:,1);

plot(v(:), j(:)*1000, 'color', colors_JV{2}, 'LineWidth', 3) 
hold on
plot(v_el(1:151), j_el(1:151)*1000, '--', 'color', colors_JV{2}, 'LineWidth', 3)
plot(v(:), zeros(1,length(v)), 'black', 'LineWidth', 1)
hold off

set(gca, 'FontSize', 25)
xlim([0, 1.2])
ylim([0,27])
xlabel('Voltage (V)', 'FontSize', 30)
ylabel('Current Density (mAcm^{-2})', 'FontSize', 30)
legend({'  Mobile Ions', '  No Mobile Ions'}, 'Location', 'southwest', 'FontSize', 30)
ax1 = gca;
%% Break down contributions to the current
%Columns in J_values are J_gen, J_rad, J_srh, J_vsr (left), J_vsr (right) and J_ext

num_values = length(CV_solutions_ion{1}.t);
J_values = zeros(num_values,7,num_devices);
J_values_el = zeros(num_values,7,num_devices);
Voc_ion = zeros(num_devices,2);
Voc_el = zeros(num_devices,2);
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

    Voc_ion(k,1) = CVstats(CVsol).Voc_f;
    Voc_ion(k,2) = interp1(v(1:ceil(num_values/2)), J_values(1:ceil(num_values/2),2,k), Voc_ion(k,1));

end    

for k=1:num_devices
    CVsol = CV_solutions_el{k};
    loss_currents = dfana.calcr(CVsol,'sub');
    x = CVsol.par.x_sub;
    num_points = length(x);
    gxt = dfana.calcg(CVsol);
    J = dfana.calcJ(CVsol);
    j_surf_rec = dfana.calcj_surf_rec(CVsol);

    J_values_el(:,1,k) = e*trapz(x, gxt(1,:))';
    J_values_el(:,2,k) = e*trapz(x, loss_currents.btb, 2)';
    J_values_el(:,3,k) = e*trapz(x, loss_currents.srh, 2)';
    J_values_el(:,4,k) = e*trapz(x(1:num_points), loss_currents.vsr(:,1:num_points), 2)';
    J_values_el(:,5,k) = e*trapz(x(num_points+1:end), loss_currents.vsr(:,num_points+1:end), 2)';
    J_values_el(:,6,k) = e*(j_surf_rec.tot);
    J_values_el(:,7,k) = J.tot(:,1);

    Voc_el(k,1) = CVstats(CVsol).Voc_f;
    Voc_el(k,2) = interp1(v(1:ceil(num_values/2)), J_values_el(1:ceil(num_values/2),2,k), Voc_el(k,1));
    
end   

%% Plot contributons to the current
%J_rad not corrected for EL - see EL_Measurements
figure(333)
num=1;
line_colour = {[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560]...
                [0 0.4470 0.7410], [0.3010 0.7450 0.9330], 'black',[0.4660 0.6740 0.1880]};
V = dfana.calcVapp(CV_solutions_ion{1});
for n = 1:7
    if n ==2 
    plot(V(:), J_values(:,n,num)*100, 'color', line_colour{n})
    else
    plot(V(:), J_values(:,n,num), 'color', line_colour{n})
    end
    hold on
end
plot(V(1:num_values), zeros(1,num_values), 'black', 'LineWidth', 1)
hold off
xlim([0, 1.2])
xlabel('Voltage (V)')
ylim([-0.03, 0.01])
ylabel('Current Density (Acm^{-2})')
legend({'J_{gen}', 'J_{rad}x100', 'J_{SRH}', 'J_{interface (left)}', 'J_{interface (right)}', '','J_{ext}'}, 'Location', 'bestoutside')

%% Losses at SC, ions vs no ions
x = categorical({'With Mobile Ions', 'Without Mobile Ions'});
y = [-J_values(31,2,1) -J_values(31,3,1) -J_values(31,4,1) -J_values(31,5,1); ...
    -J_values_el(31,2,1) -J_values_el(31,3,1) -J_values_el(31,4,1) -J_values_el(31,5,1)];

bar_colours = {[0.9290 0.6940 0.1250],[0 0.4470 0.7410],[0.4660 0.6740 0.1880],[0.3010 0.7450 0.9330]};
figure('Name', 'RecombinationSC_ElvsIon','Position', [100 100 1250 2000])
box on
b = bar(x, y, 'stacked','FaceColor','flat');

for k = 1:4
    b(k).CData = bar_colours{k};
end

set(gca, 'FontSize', 30)
ylabel('Recombination Current Density (Acm^{-2})', 'FontSize', 30)
legend({'', '  J_{bulk}', '  J_{surface}', ''}, 'Location', 'northeast', 'FontSize', 30)
%% Plot PLQY results
%One plot for effect of changing mobility 
%Another plot for effect of changing energetics
figure('Name', 'PLQYPlot Mobility', 'Position', [100 100 1250 2000])
mob_order = [2,3,1];
for i = 1:3
    j = mob_order(i);
    semilogy(V(1:151), 100*(J_values(1:151,2,i))./J_values(1:151,1,i), 'color', colors_JV{j}, 'LineWidth', 3) 
    hold on 
    semilogy(Voc_ion(i,1), 100*Voc_ion(i,2)/J_values(70,1,i), 'ko', 'MarkerSize', 15, 'LineWidth', 2)
    semilogy(V(1:151), 100*(J_values_el(1:151,2,i))./J_values_el(1:151,1,i), '--', 'color', colors_JV{j}, 'LineWidth', 3) 
    semilogy(Voc_el(i,1), 100*Voc_el(i,2)/J_values_el(70,1,i), 'ko', 'MarkerSize', 15, 'LineWidth', 2)
end
hold off

set(gca, 'FontSize', 25)
xlim([0, 1.2])
ylim([1e-5, 0.3])
legend({'  5 x 10^{-4}', '', '','','  5 x 10^{-5}', '', '', '', '  5 x 10^{-6}', '', '', ''}, 'Location', 'southeast', 'FontSize', 30)
xlabel('Voltage (V)', 'FontSize', 30)
ylabel('PLQY (%)', 'FontSize', 30)
title(legend, 'ETM Mobility (cm^{2}/Vs)', 'Fontsize', 30)
yyaxis right
kT = 0.0257;
J_gen = -J_values(70,1,1);
ni = CV_solutions_ion{1}.par.ni(3);
B = CV_solutions_ion{1}.par.B(3);
w = CV_solutions_ion{1}.par.d(3);
J0 = -e*w*B*ni^2;
ylim([kT*log((1e-7*J_gen)/J0), kT*log((0.003*J_gen)/J0)])
ylabel('QFLS (eV)', 'FontSize', 30)
yticks([0.90, 0.95, 1.00, 1.05, 1.10])
yticklabels({'0.90', '0.95', '1.00', '1.05', '1.10'})
ax2=gca;
ax2.YAxis(1).Color = 'k';
ax2.YAxis(2).Color = 'k';

figure('Name', 'PLQYPlot Energetics', 'Position', [100 100 1250 2000])
energetics_order_devices = [4,1,5];
for j = 1:3
    i = energetics_order_devices(j);
    semilogy(V(1:151), 100*(J_values(1:151,2,i))./J_values(1:151,1,i), 'color', colors_JV{j}, 'LineWidth', 3) 
    hold on 
    semilogy(Voc_ion(i,1), 100*Voc_ion(i,2)/J_values(70,1,i), 'ko', 'MarkerSize', 15, 'LineWidth', 2)
    semilogy(V(1:151), 100*(J_values_el(1:151,2,i))./J_values_el(1:151,1,i), '--', 'color', colors_JV{j}, 'LineWidth', 3) 
    semilogy(Voc_el(i,1), 100*Voc_el(i,2)/J_values_el(70,1,i), 'ko', 'MarkerSize', 15, 'LineWidth', 2)
end
hold off

set(gca, 'FontSize', 25)
xlim([0, 1.2])
ylim([1e-5, 0.3])
legend({'  4.15', '', '', '', '  3.95', '', '', '', '  3.75', '', '', ''}, 'Location', 'southeast', 'FontSize', 30)
xlabel('Voltage (V)', 'FontSize', 30)
ylabel('PLQY (%)', 'FontSize', 30)
title(legend, '|ETM LUMO| (eV)', 'Fontsize', 30)
yyaxis right
ylim([kT*log((1e-7*J_gen)/J0), kT*log((0.003*J_gen)/J0)])
ylabel('QFLS (eV)', 'FontSize', 30)
yticks([0.90, 0.95, 1.00, 1.05, 1.10])
yticklabels({'0.90', '0.95', '1.00', '1.05', '1.10'})
ax3=gca;
ax3.YAxis(1).Color = 'k';
ax3.YAxis(2).Color = 'k';

%% Save Plots at 300 dpi
save = 1;
fig_num = 3;

if save == 1 && fig_num == 1
    exportgraphics(ax1, ...
    'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Weidong_ETL\Paper\v2\JV.png', ...
    'Resolution', 300)
elseif save == 1 && fig_num == 2
    exportgraphics(ax2, ...
    'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Weidong_ETL\Paper\v2\PLQY_mob.png', ...
    'Resolution', 300)
elseif save == 1 && fig_num == 3
    exportgraphics(ax3, ...
    'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Weidong_ETL\Paper\v2\PLQY_energetics.png', ...
    'Resolution', 300)
end 










