%Code to model the results of Weidong's experiment varying the choice of
%ETL material used in his FaCs devices. Just without the FaCs and the worng
%HTL...

%% Read in data files 

par_kloc6 = pc('Input_files/pedotpss_mapi_kloc6.csv');
par_icba = pc('Input_files/pedotpss_mapi_icba.csv');
par_pcbm = pc('Input_files/pedotpss_mapi_pcbm.csv');

devices = {par_kloc6, par_pcbm, par_icba};

%% Find SS solutions at 1 sun intensity (for SC QFLS)
eqm_solutions_dark = cell(1,3);
eqm_solutions_light = cell(1,3);
for i = 1:3
    eqm = equilibrate(devices{i});
    eqm_solutions_dark{i} = eqm;
    eqm_solutions_light{i} = changeLight(eqm.ion,1,0);
end

%% Perform CV scans
CV_solutions_el = cell(1,3);
CV_solutions_ion = cell(1,3);
for j = 1:3
    sol_el = eqm_solutions_dark{j}.el;
    sol_ion = eqm_solutions_dark{j}.ion;
    CV_solutions_el{j} = doCV(sol_el, 1, -0.3, 1.3, -0.3, 1e-3, 1, 321);
    CV_solutions_ion{j} = doCV(sol_ion, 1, -0.3, 1.3, -0.3, 1e-3, 1, 321);
end

%% Plot JVs
figure(1)
for m=1:3
    v = dfana.calcVapp(CV_solutions_ion{m});
    j = dfana.calcJ(CV_solutions_ion{m}).tot(:,1);
    plot(v(:), j(:))
    hold on
end
plot(v(:), zeros(1,length(v)), 'black', 'LineWidth', 1)
hold off
legend({'Kloc-6', 'PCBM', 'ICBA',''}, 'Location', 'northwest')
xlim([0, 1.3])
ylim([-0.03, 0.01])
xlabel('Voltage(V)')
ylabel('Current Density (Acm^{-2})')

%% Break down contributions to the current
%Columns in J_values are J_gen, J_rad, J_srh, J_vsr and J_ext
CV_solutions = CV_solutions_ion;

num_values = length(CV_solutions{1}.t);
J_values = zeros(num_values, 6,3);
e = -CV_solutions{1}.par.e;

for k=1:3
    CVsol = CV_solutions{k};
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
%% Plot contributons to the current
figure(2)
line_colour = {[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560],...
                [0 0.4470 0.7410], [0.3010 0.7450 0.9330], [0.4660 0.6740 0.1880]};
V = dfana.calcVapp(CV_solutions{1});
for n = 1:6
    plot(V(:), J_values(:,n,1), 'color', line_colour{n})
    hold on
end
plot(V(1:num_values), zeros(1,num_values), 'black', 'LineWidth', 1)
hold off
xlim([0, max(V)])
xlabel('Voltage (V)')
ylim([-0.025, 0.01])
ylabel('Current Density (Acm^{-2})')
legend({'J_{gen}', 'J_{rad}', 'J_{SRH}', 'J_{VSR}', 'J_{surf}', 'J_{ext}'}, 'Location', 'bestoutside')

%% Plot 'PL' results
figure(3)
for i = 1:3
    plot(dfana.calcVapp(CV_solutions{i}), -J_values(:,2,i)) 
    hold on
end
xlim([0, 1.3])
xlabel('Voltage (V)')
ylabel('e\phi_{PL}(Acm^{-2})')
ylim([0, 0.5e-3])
legend({'Kloc-6', 'PCBM', 'ICBA'}, 'Location', 'northwest')

%% Plot rad/non-rad ratio
figure(4)
for n = 1:3
    semilogy(dfana.calcVapp(CV_solutions{n}), J_values(:,2,n)./(J_values(:,3,n)+J_values(:,4,n)+J_values(:,5,n)))
    hold on
end
xlim([0, 1.3])
xlabel('Voltage (V)')
ylabel('J_{rad}/J_{non-rad}')
%ylim([0, 2.5e-3])
legend({'Kloc-6', 'PCBM', 'ICBA'}, 'Location', 'southwest')