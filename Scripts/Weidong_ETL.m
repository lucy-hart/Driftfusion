%Code to model the results of Weidong's experiment varying the choice of
%ETL material used in his FaCs devices.

%% Read in data files 
num_devices = 3;

par_kloc6 = pc('Input_files/PTAA_MAPI_Kloc6_v4.csv');
par_pcbm = pc('Input_files/PTAA_MAPI_PCBM_v4.csv');
par_icba = pc('Input_files/PTAA_MAPI_ICBA_v4.csv');
par_iph = pc('Input_files/PTAA_MAPI_IPH_v4.csv');

if num_devices == 4
    devices = {par_kloc6, par_pcbm, par_icba, par_iph};
elseif num_devices == 3
    devices = {par_kloc6, par_pcbm, par_icba};
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
    suns = 1.15;
    CV_solutions_el{j} = doCV(eqm_solutions_dark{j}.el, suns, -0.3, 1.2, -0.3, 10e-3, 1, 301);
    CV_solutions_ion{j} = doCV(eqm_solutions_dark{j}.ion, suns, -0.3, 1.2, -0.3, 10e-3, 1, 301);
end

%% Plot JVs
figure('Name', 'JVPlot', 'Position', [100 100 1250 1250])
colors_JV = {[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],[0 0.4470 0.7410],[0.9290 0.6940 0.1250]};

for m=1:num_devices
    v = dfana.calcVapp(CV_solutions_ion{m});
    v_el = dfana.calcVapp(CV_solutions_el{m});
    j = -dfana.calcJ(CV_solutions_ion{m}).tot(:,1);
    j_el = -dfana.calcJ(CV_solutions_el{m}).tot(:,1);
    plot(v(:), j(:)*1000, 'color', colors_JV{m}, 'LineWidth', 3) 
    hold on
    plot(v_el(1:151), j_el(1:151)*1000, '--', 'color', colors_JV{m}, 'LineWidth', 3)
    hold on
end
plot(v(:), zeros(1,length(v)), 'black', 'LineWidth', 1)
hold off

set(gca, 'FontSize', 20)
xlim([0, 1.2])
ylim([0,27])
legend({'ETM 1','', 'ETM 2','', 'ETM 3','',''}, 'Location', 'southwest', 'FontSize', 30)
xlabel('Voltage(V)', 'FontSize', 30)
ylabel('Current Density (mAcm^{-2})', 'FontSize', 30)
ax1 = gca;
%% Break down contributions to the current
%Columns in J_values are J_gen, J_rad, J_srh, J_vsr and J_ext
%Don't take btb from back layer as unlikely to escape the device, parasitic
%absorption of metal back contact

num_values = length(CV_solutions_ion{1}.t);
num_stop = sum(CV_solutions_ion{1}.par.layer_points(1:3));
J_values = zeros(num_values,7,num_devices);
J_values_el = zeros(num_values,2,num_devices);
e = -CV_solutions_ion{1}.par.e;

for k=1:num_devices
    CVsol = CV_solutions_ion{k};
    loss_currents = dfana.calcr(CVsol,'sub');
    x = CVsol.par.x_sub;
    gxt = dfana.calcg(CVsol);
    J = dfana.calcJ(CVsol);
    j_surf_rec = dfana.calcj_surf_rec(CVsol);

    J_values(:,1,k) = e*trapz(x, gxt(1,:))';
    J_values(:,2,k) = e*trapz(x(1:num_stop), loss_currents.btb(:,1:num_stop), 2)';
    J_values(:,3,k) = e*trapz(x, loss_currents.srh, 2)';
    J_values(:,4,k) = e*trapz(x, loss_currents.vsr, 2)';
    J_values(:,5,k) = e*(j_surf_rec.tot);
    J_values(:,6,k) = J.tot(:,1);
    
end    

for k=1:num_devices
    CVsol = CV_solutions_el{k};
    loss_currents = dfana.calcr(CVsol,'sub');
    x = CVsol.par.x_sub;
    gxt = dfana.calcg(CVsol);

    J_values_el(:,1,k) = e*trapz(x, gxt(1,:))';
    J_values_el(:,2,k) = e*trapz(x(1:num_stop), loss_currents.btb(:,1:num_stop), 2)';
    
end   

%% Plot contributons to the current
%J_rad not corrected for EL - see EL_Measurements
figure(333)
num=3;
line_colour = {[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560]...
                [0 0.4470 0.7410], [0.3010 0.7450 0.9330], [0.4660 0.6740 0.1880]};
V = dfana.calcVapp(CV_solutions_ion{1});
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

%% Plot PLQY results
figure('Name', 'PLQYPlot', 'Position', [100 100 1250 2000])
for i = 1:num_devices
    semilogy(V(1:151), 100*(J_values(1:151,2,i))./J_values(1:151,1,i), 'color', colors_JV{i}, 'LineWidth', 3) 
    hold on 
    semilogy(v_el(1:151), 100*(J_values_el(1:151,2,i))./J_values_el(1:151,1,i), '--', 'color', colors_JV{i}, 'LineWidth', 3)     
    hold on
end
hold off

set(gca, 'FontSize', 20)
xlim([0, 1.2])
ylim([1e-5, 0.3])
legend({'ETM 1','', 'ETM 2','', 'ETM 3',''}, 'Location', 'southeast', 'FontSize', 30)
xlabel('Voltage (V)', 'FontSize', 30)
ylabel('PLQY (%)', 'FontSize', 30)
ax2=gca;

%% Save Plots at 300 dpi
save = 1;
fig_num = 2;

if save == 1 && fig_num == 1
    exportgraphics(ax1, ...
    'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Weidong_ETL\Paper\JV.png', ...
    'Resolution', 300)
elseif save == 1 && fig_num == 2
    exportgraphics(ax2, ...
    'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Weidong_ETL\Paper\PLQY.png', ...
    'Resolution', 300)
end 

%% Plot Jnonrad/Jrad results
figure('Name', 'Jnr/Jrad', 'Position', [100 100 1250 2000])
for i = 1:num_devices
    semilogy(V(1:151), (sum(J_values(1:151,3:5,i),2)./(J_values(1:151,2,i)))./(sum(J_values(31,3:5,i),2)./(J_values(31,2,i))) ...
        , 'color', colors_JV{i}, 'LineWidth', 3) 
    hold on 
end
hold off

set(gca, 'FontSize', 20)
xlim([0, 1.2])
ylim([0.3, 3])
legend({'ETM 1','ETM 2','ETM 3'}, 'Location', 'southeast', 'FontSize', 30)
xlabel('Voltage (V)', 'FontSize', 30)
ylabel('J_{non rad}/J_{rad} (norm)', 'FontSize', 30)
ax2=gca;








