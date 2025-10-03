par = pc('Input_files/SAM_MAFACsPbIBr_C60.csv');
par.vsr_mode = 0;
par.Rs = 1e6;
par = refresh_device(par);
eqm = equilibrate(par,0,0);

%%
eqm.ion.par.laser_lambda2 = 605;
% eqm.ion.par.g2_fun_type = 'square';
% eqm.ion.par.g2_fun_arg = [0 power t_stop duty];
eqm.ion.par.g2_fun_type = 'constant';
eqm.ion.par.g2_fun_arg = [power];
eqm.ion.par.int2 = [1];
eqm.ion.par.g1_fun_arg = [0];
eqm.ion.par.tmax = 1;
eqm.ion.par.tmesh_type = 'linear';
eqm.ion.par.side = 'right';
eqm.ion.par = refresh_device(eqm.ion.par);
sol_light = df(eqm.ion);
%%
t_stop = 1e-6;
%As a percentage
duty = 0.1;
power = 10;
sol_light.par.tmax = 0.95*t_stop;
sol_light.par.tpoints = 5e3;
sol_light.par.PPP = 1;
sol_light.par.PPP_args = [0 1e4 t_stop duty];
sol_light.par = refresh_device(sol_light.par);

pulse_sol_kinetic = df(sol_light);

t_kinetic = pulse_sol_kinetic.t;
PL_kinetic = dfana.calcPLt(pulse_sol_kinetic);

%%
J = dfana.calcJ(pulse_sol_kinetic);
t_start_mean = 0.3*(duty/100)*t_stop;
t_end_mean = 0.4*(duty/100)*t_stop;
[~,arg_start_mean] = min(abs(pulse_sol_kinetic.t - t_start_mean));
[~,arg_end_mean] = min(abs(pulse_sol_kinetic.t - t_end_mean));

J_lighton = mean(J.tot(arg_start_mean:arg_end_mean,1));


%%
figure('Name', 'PPP Current')

plot(pulse_sol_kinetic.t(arg_start_mean:end), J.tot(arg_start_mean:end,1)-J_lighton)
xlabel('Time (s)')
ylabel('PPP Photocurrent (A cm-2)')

%%
figure('Name', 'Total Current')

plot(pulse_sol_kinetic.t, J.tot(:,42))
xlabel('Time (s)')
ylabel('Total Photocurrent (A cm-2)')
