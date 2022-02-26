%Code to model the results of Weidong's experiment varying the choice of
%ETL material used in his FaCs devices. Just without the FaCs..

%% Read in data files 

par_kloc6 = pc('Input_files/PTAA_MAPI_Kloc6_v2.csv');
par_pcbm = pc('Input_files/PTAA_MAPI_PCBM_v2.csv');
par_icba = pc('Input_files/PTAA_MAPI_ICBA_v2.csv');

devices = {par_kloc6, par_pcbm, par_icba};

%% Find eqm solutions 
eqm_solutions_dark = cell(1,3);
for i = 1:3
    eqm_solutions_dark{i} = equilibrate(devices{i});
end

%% Perform CV scans
CV_solutions_el = cell(1,3);
CV_solutions_ion = cell(1,3);
for j = 1:3
    CV_solutions_el{j} = doCV(eqm_solutions_dark{j}.el, 1.2, -0.3, 1.3, -0.3, 1e-3, 1, 321);
    CV_solutions_ion{j} = doCV(eqm_solutions_dark{j}.ion, 1.2, -0.3, 1.3, -0.3, 1e-3, 1, 321);
end

%% Plot JVs
figure(1)
colors_JV = {[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],[0 0.4470 0.7410]};
for m=1:3
    v = dfana.calcVapp(CV_solutions_ion{m});
    v_el = dfana.calcVapp(CV_solutions_el{m});
    j = -dfana.calcJ(CV_solutions_ion{m}).tot(:,1);
    j_el = -dfana.calcJ(CV_solutions_el{m}).tot(:,1);
    plot(v(:), j(:), 'color', colors_JV{m}) 
    hold on
    plot(v_el(1:161), j_el(1:161), '-.', 'color', colors_JV{m})
    hold on
end
plot(v(:), zeros(1,length(v)), 'black', 'LineWidth', 1)
hold off
legend({'Kloc-6','', 'PCBM','', 'ICBA','',''}, 'Location', 'southwest')
xlim([0, 1.3])
ylim([0, 0.03])
xlabel('Voltage(V)')
ylabel('Current Density (Acm^{-2})')

%% Break down contributions to the current
%Columns in J_values are J_gen, J_rad, J_srh, J_vsr and J_ext
%Don't take btb from back layer as unlikely to escape the device, parasitic
%absorption of metal back contact

num_values = length(CV_solutions_ion{1}.t);
num_stop = sum(CV_solutions_ion{1}.par.layer_points(1:3));
J_values = zeros(num_values,7,3);
J_values_el = zeros(num_values,2,3);
e = -CV_solutions_ion{1}.par.e;

for k=1:3
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

for k=1:3
    CVsol = CV_solutions_el{k};
    loss_currents = dfana.calcr(CVsol,'sub');
    x = CVsol.par.x_sub;
    gxt = dfana.calcg(CVsol);

    J_values_el(:,1,k) = e*trapz(x, gxt(1,:))';
    J_values_el(:,2,k) = e*trapz(x(1:num_stop), loss_currents.btb(:,1:num_stop), 2)';
    
end   
%% Plot contributons to the current
%J_rad not corrected for EL - see EL_Measurements
figure(3)
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
xlim([0, 1.3])
xlabel('Voltage (V)')
ylim([-0.03, 0.01])
ylabel('Current Density (Acm^{-2})')
legend({'J_{gen}', 'J_{rad}x100', 'J_{SRH}', 'J_{interface}', 'J_{contact}','J_{ext}'}, 'Location', 'bestoutside')

%% Plot PLQY results
figure(4)
for i = 1:3
    semilogy(V(1:161), 100*(J_values(1:161,2,i))./J_values(1:161,1,i), 'color', colors_JV{i}) 
    hold on 
    semilogy(V(1:161), 100*(J_values_el(1:161,2,i))./J_values_el(1:161,1,i), '-.', 'color', colors_JV{i})     
    hold on
end
hold off
xlim([0, 1.3])
xlabel('Voltage (V)')
ylabel('PLQY (%)')
ylim([1e-4, 0.3])
legend({'Kloc-6','', 'PCBM','', 'ICBA',''}, 'Location', 'northwest')





