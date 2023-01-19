function Plot_Current_Contributions(CVsol)
%Script to plot the currents from radiative and non-radiative (split into
%SRH, surface and VSR) losses, alongside the current measured in a JV sweep and the
%generation current. 

%% Break down contributions to the current
%Columns in J_values are J_gen, J_rad, J_srh, J_vsr, J_shunt and J_ext
num_values = length(CVsol.t);
J_values = zeros(num_values, 6);
e = -CVsol.par.e;
loss_currents = dfana.calcr(CVsol,'sub');
x = CVsol.par.x_sub;
gxt = dfana.calcg(CVsol);
J = dfana.calcJ(CVsol);
j_surf_rec = dfana.calcj_surf_rec(CVsol);

%forward sweep
J_values(:,1) = e*trapz(x, gxt(1,:))';
J_values(:,2) = 100*e*trapz(x, loss_currents.btb, 2)';
J_values(:,3) = e*trapz(x, loss_currents.srh, 2)';
J_values(:,4) = e*trapz(x, loss_currents.vsr, 2)';
J_values(:,5) = e*(j_surf_rec.tot);
J_values(:,6) = J.tot(:,1);
    
%% Plot contributons to the current
figure('Name', 'Current Contributions', 'Position', [50 50 1000 1000])
line_colour = {[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560],...
                [0 0.4470 0.7410], [0.3010 0.7450 0.9330], [0.4660 0.6740 0.1880]};
V = dfana.calcVapp(CVsol);
for n = 1:6
    plot(V(:), J_values(:,n), 'color', line_colour{n})
    hold on
end
plot(V(1:num_values), zeros(1,num_values), 'black', 'LineWidth', 1)
hold off
V1 = CVsol.par.V_fun_arg(2);
V2 = CVsol.par.V_fun_arg(3);
xlim([0, max([V1, V2])])
xlabel('Voltage (V)')
ylim([J_values(1,1)*1.1, 0.01])
ylabel('Current Density (Acm^{-2})')
legend({'J_{gen}', 'J_{rad}x100', 'J_{SRH}', 'J_{surface}', '', 'J_{ext}'}, 'Location', 'bestoutside')

end