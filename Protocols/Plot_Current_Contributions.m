function Plot_Current_Contributions(CVsol)
%Script to plot the currents from radiative and non-radiative (split into SRH
%and VSR) losses, alongside the current measured in a JV sweep and the
%generation current. 

%% Break down contributions to the current
%Columns in J_values are J_gen, J_rad, J_srh, J_vsr and J_ext
num_values = length(CVsol.t);
J_values = zeros(num_values,5);
e = CVsol.par.e;
loss_currents = dfana.calcr(CVsol,'whole');

%forward sweep
J_values(:,1) = -e*trapz(CVsol.par.x_sub, dfana.calcg(CVsol).').';
J_values(:,2) = -e*trapz(CVsol.par.x_sub, loss_currents.btb(:,1:end-1).').';
J_values(:,3) = -e*trapz(CVsol.par.x_sub, loss_currents.srh(:,1:end-1).').';
J_values(:,4) = -e*trapz(CVsol.par.x_sub, loss_currents.vsr(:,1:end-1).').';
J_values(:,5) = dfana.calcJ(CVsol).tot(:,1);
    
%% Plot contributons to the current
figure(2)
line_colour = {[0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0 0.4470 0.7410], [0.3010 0.7450 0.9330], [0.4660 0.6740 0.1880]};
V = dfana.calcVapp(CVsol);
hold on
for n = 1:5
    plot(V(:), J_values(:,n), 'color', line_colour{n}) 
end
plot(V(1:num_values), zeros(1,num_values), 'black', 'LineWidth', 1)
hold off
xlim([0, CVsol.par.V_fun_arg(2)])
xlabel('Voltage (V)')
ylim([-0.025, 0.025])
ylabel('Current Density (Acm^{-2})')
legend({'J_{gen}', 'J_{rad}', 'J_{SRH}', 'J_{VSR}', 'J_{ext}', ''}, 'Location', 'bestoutside')