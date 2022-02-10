%Simulate a CE at SC experiment

%% Load in device (assuming I've run Weidong_ETL first)
eqm_CE = devices{2};
eqm_CE = equilibrate(eqm_CE).ion;

%% Make variables
num_start = sum(eqm_CE.par.layer_points(1:2))+1;
num_stop = num_start + eqm_CE.par.layer_points(3)-1;
x = eqm_CE.par.x_sub;

light_intensities = logspace(-1,1,10);
Q_SC_meas = zeros(1,length(light_intensities));
transient_currents = cell(1,length(light_intensities));
Q_SC_act = zeros(1,length(light_intensities));

%%
j=1;
for i=light_intensities   
    SC_sol = jumptoV(eqm_CE, 0, 1e7, 1, i, 1, 0);
    n_SC = trapz(x(num_start:num_stop), SC_sol.u(end, num_start:num_stop, 2), 2)*0.045;
    p_SC = trapz(x(num_start:num_stop), SC_sol.u(end, num_start:num_stop, 3), 2)*0.045;
    Q_SC_act(j) = -e*(n_SC*p_SC)^0.5;
    SC_sol.par.int1=0;
    SC_sol.par.t0 = 1e-15;
    SC_sol.par.tmax = 2e-5;
    SC_sol.par.tpoints = 500;
    SC_sol.par.tmesh_type = 'log10';
    CE_transient = df(SC_sol);
    transient_currents{j} = dfana.calcJ(CE_transient);
    J_ex = transient_currents{j}.tot(:,1);
    Q_SC_meas(j) = trapz(CE_transient.t, J_ex)*0.045;
    j=j+1;
end

%% Plot

figure(7)
loglog(light_intensities, -Q_SC_meas, 'o', light_intensities, Q_SC_act, 'x')

figure(8)
for k = 1:length(light_intensities)
    plot(CE_transient.t, -transient_currents{k}.tot(:,1))
    hold on
end
hold off 

xlim([0, 2e-7])