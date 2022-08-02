%ESA code
%Investigate effect of ETM energetics and ion motion on device performance
%Run JV scans for all 3 ETMs and plot the contributons to the current

%% Read in data files 
par_ETM1 = pc('Input_files/PTAA_MAPI_NegOffset.csv');
par_ETM2 = pc('Input_files/PTAA_MAPI_NoOffset.csv');
par_ETM3 = pc('Input_files/PTAA_MAPI_PosOffset.csv');

num_devices = 3;
devices = {par_ETM1, par_ETM2, par_ETM3};

%% Find eqm solutions 
eqm_solutions_dark = cell(1,num_devices);
for i = 1:num_devices
    %Supress error message about vsr for ETM1 (the error is 5.2% and the warning limit threshold in 5.0%)
    devices{i}.RelTol_vsr = 0.06;
    devices{i} = refresh_device(devices{i});
    eqm_solutions_dark{i} = equilibrate(devices{i});
end

%% Perform JV scans
JV_solutions_el = cell(1,num_devices);
JV_solutions_ion = cell(1,num_devices);
suns = 1.15;
for j = 1:num_devices
    JV_solutions_el{j} = doCV(eqm_solutions_dark{j}.el, suns, -0.3, 1.2, -0.3, 10e-3, 1, 301);
    JV_solutions_ion{j} = doCV(eqm_solutions_dark{j}.ion, suns, -0.3, 1.2, -0.3, 10e-3, 1, 301);
end

%% Plot JVs
figure('Name', 'JVPlot', 'Position', [100 100 1250 1250])
colors_JV = {[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],[0 0.4470 0.7410],[0.9290 0.6940 0.1250]};

for m=1:num_devices
    v = dfana.calcVapp(JV_solutions_ion{m});
    v_el = dfana.calcVapp(JV_solutions_el{m});
    j = dfana.calcJ(JV_solutions_ion{m}).tot(:,1);
    j_el = dfana.calcJ(JV_solutions_el{m}).tot(:,1);
    plot(v(:), j(:)*1000, 'color', colors_JV{m}, 'LineWidth', 3) 
    hold on
    plot(v_el(1:151), j_el(1:151)*1000, '--', 'color', colors_JV{m}, 'LineWidth', 3)
    hold on
end
plot(v(:), zeros(1,length(v)), 'black', 'LineWidth', 1)
hold off

set(gca, 'FontSize', 25)
xlim([-0.2, 1.2])
ylim([-27,8])
legend({'ETM 1','', 'ETM 2','', 'ETM 3','',''}, 'Location', 'northwest', 'FontSize', 30)
xlabel('Voltage(V)', 'FontSize', 30)
ylabel('Current Density (mAcm^{-2})', 'FontSize', 30)
ax1 = gcf;
%% Break down contributions to the current
%Columns in J_values are J_gen, J_rad, J_srh, J_vsr and J_ext

num_values = length(JV_solutions_ion{1}.t);
num_stop = sum(JV_solutions_ion{1}.par.layer_points(1:3));
J_values = zeros(num_values,7,num_devices);
J_values_el = zeros(num_values,2,num_devices);
e = -JV_solutions_ion{1}.par.e;

for k=1:num_devices
    CVsol = JV_solutions_ion{k};
    loss_currents = dfana.calcr(CVsol,'sub');
    x = CVsol.par.x_sub;
    gxt = dfana.calcg(CVsol);
    J = dfana.calcJ(CVsol);
    j_surf_rec = dfana.calcj_surf_rec(CVsol);

    J_values(:,1,k) = e*trapz(x, gxt(1,:))';
    J_values(:,2,k) = e*trapz(x, loss_currents.btb, 2)';
    J_values(:,3,k) = e*trapz(x, loss_currents.srh, 2)';
    J_values(:,4,k) = e*trapz(x, loss_currents.vsr, 2)';
    J_values(:,5,k) = e*(j_surf_rec.tot);
    J_values(:,6,k) = J.tot(:,1);
    
end    

for k=1:num_devices
    CVsol = JV_solutions_el{k};
    loss_currents = dfana.calcr(CVsol,'sub');
    x = CVsol.par.x_sub;
    gxt = dfana.calcg(CVsol);

    J_values_el(:,1,k) = e*trapz(x, gxt(1,:))';
    J_values_el(:,2,k) = e*trapz(x, loss_currents.btb, 2)';
    
end   

%% Plot contributons to the current
%J_rad not corrected for EL - see EL_Measurements
figure(333)
num=2;
line_colour = {[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560]...
                [0 0.4470 0.7410], [0.3010 0.7450 0.9330], [0.4660 0.6740 0.1880]};
V = dfana.calcVapp(JV_solutions_ion{1});
for n = 1:6
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
legend({'J_{gen}', 'J_{rad}x100', 'J_{SRH}', 'J_{interface}', 'J_{contact}','J_{ext}'}, 'Location', 'bestoutside')
ax2 = gca;
%% Save Plots at 300 dpi
save = 1;
fig_num = 1;
PS = 1;

if save == 1 && fig_num == 1
    exportgraphics(ax1, ...
    'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\ESA\Figures\JV.png', ...
    'Resolution', 300)
elseif save == 1 && fig_num == 2
    filename = ['C:\Users\ljh3218\OneDrive - Imperial College London\PhD\ESA\Figures\CurrentContributions_PS' num2str(PS_distplot) '.png'];
    exportgraphics(ax2, filename, 'Resolution', 300)
end 







