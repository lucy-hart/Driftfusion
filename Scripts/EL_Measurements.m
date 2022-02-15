%Program to do the EL measurements from Weidong's dataset.
%Run this after Weidong_ETL

%% Do dark JV
dark_CV_solutions = cell(1,3);
for i = 1:3
    dark_CV_solutions{i} = doCV(eqm_solutions_dark{i}.ion, 0, -0.25, 1.25, -0.25, 1e-3, 1, 321);
end
%% Plot dark JVs
figure(11)
for m=1:3
    v_dark = dfana.calcVapp(dark_CV_solutions{m});
    j_dark = dfana.calcJ(dark_CV_solutions{m}).tot(:,1);
    plot(v_dark(:), j_dark(:))
    hold on
end
plot(v_dark(:), zeros(1,length(v_dark)), 'black', 'LineWidth', 1)
hold off
legend({'Kloc-6', 'PCBM', 'ICBA',''}, 'Location', 'southwest')
xlim([0, 1.3])
ylim([0, 0.025])
xlabel('Voltage(V)')
ylabel('Current Density (Acm^{-2})')

%% Break down contributions to the current
%Columns in J_values_dark are J_EL, J_srh, J_vsr, J_ext and J_surf
num_values = length(dark_CV_solutions{1}.t);
J_values_dark = zeros(num_values, 5,3);
e = dark_CV_solutions{1}.par.e;

for k=1:3
    CVsol = dark_CV_solutions{k};
    loss_currents = dfana.calcr(CVsol,'sub');
    x = CVsol.par.x_sub;
    J = dfana.calcJ(CVsol);
    j_surf_rec = dfana.calcj_surf_rec(CVsol);

    J_values_dark(:,1,k) = e*trapz(x, loss_currents.btb, 2)';
    J_values_dark(:,2,k) = e*trapz(x, loss_currents.srh, 2)';
    J_values_dark(:,3,k) = e*trapz(x, loss_currents.vsr, 2)';
    J_values_dark(:,4,k) = J.tot(:,1);
    J_values_dark(:,5,k) = e*(j_surf_rec.tot);
end  

%% Plot contributons to the current corrected for EL
figure(12)
num = 2;

V = dfana.calcVapp(dark_CV_solutions{1});
plot(V(:), J_values(:,1,num), 'color', [0.8500 0.3250 0.0980])
hold on
plot(V(:), (J_values(:,2,num)+J_values_dark(:,1,num))*100, 'color', [0.9290 0.6940 0.1250])
hold on
plot(V(:), -J_values_dark(:,1,num)*100, 'r:')
hold on
plot(V(:), J_values(:,3,num), 'color', [0 0.4470 0.7410])
hold on
plot(V(:), J_values(:,4,num), 'color', [0.3010 0.7450 0.9330])
hold on
plot(V(:), J_values(:,5,num), 'color', [0.4660 0.6740 0.1880])
hold on
plot(V(1:num_values), zeros(1,num_values), 'black', 'LineWidth', 1)
hold off

xlim([0, max(V)])
xlabel('Voltage (V)')
ylim([-0.025, 0.005])
ylabel('Current Density (Acm^{-2})')
legend({'J_{gen}', 'J_{rad}x100', 'J_{EL}x100', 'J_{SRH}', 'J_{VSR}', 'J_{ext}'}, 'Location', 'bestoutside')

%% Plot PLQY results
figure(13)
for i = 1:3
    semilogy(V(1:161), 100*(J_values(1:161,2,i))./J_values(1:161,1,i), 'color', colors_JV{i})                                                       
    hold on
end
xlim([0, 1.25])
xlabel('Voltage (V)')
ylabel('PLQY (%)')
ylim([1e-3, 0.3])
legend({'Kloc-6', 'PCBM', 'ICBA'}, 'Location', 'northwest')

%% Plot rad/non-rad ratio. 
figure(14)
for n = 1:3
    semilogy(dfana.calcVapp(CV_solutions{n}), 100*(J_values(:,2,n))...
               ./(J_values(:,1,n)-J_values(:,6,n)+J_values(:,2,n)))
    hold on
end
xlim([0, 1.3])
xlabel('Voltage (V)')
ylabel('J_{rad}/J_{non-rad} (%)')
%ylim([1e-3, 0.3])
legend({'Kloc-6', 'PCBM', 'ICBA'}, 'Location', 'northwest')
