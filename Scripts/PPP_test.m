par = pc('Input_files/SAM_MAFACsPbIBr_C60.csv');
par.vsr_mode = 0;
par.Rs = 0;
par.z_t = -1;
par = refresh_device(par);
eqm = equilibrate(par);

%%
sol = eqm.ion;
sol.par.laser_lambda2 = 750;
sol.par.g2_fun_type = 'constant';
sol.par.g2_fun_arg = [1];
sol.par.int2 = [0.1];
sol.par.g1_fun_arg = [0];
sol.par.tmax = 1;
sol.par.tmesh_type = 'linear';
sol.par.side = 'right';
sol.par = refresh_device(sol.par);
sol_light = df(sol);

%%
J = dfana.calcJ(sol_light);
t_start_mean = 0.5;
t_end_mean = 0.9;
[~,arg_start_mean] = min(abs(sol_light.t - t_start_mean));
[~,arg_end_mean] = min(abs(sol_light.t - t_end_mean));

J_lighton = mean(J.tot(arg_start_mean:arg_end_mean,1));

%%
sol = eqm.ion;
sol.par.g2_fun_arg = [0];
sol.par.int2 = [0];
sol.par.g1_fun_arg = [0];
sol.par.tmesh_type = 'linear';
sol.par.side = 'right';

period = 2e-8;
%As a percentage
duty = 50;
sol.par.tmax = 0.95*(period);
sol.par.tpoints = 1e3;
sol.par.PPP = 1;
sol.par.PPP_args = [0 1e5 period duty];
sol.par = refresh_device(sol.par);

sol_darkpulse = df(sol);
J_pulse_dark = dfana.calcJ(sol_darkpulse);
    
%%
sol = eqm.ion;
sol.par.laser_lambda2 = 750;
sol.par.g2_fun_type = 'square';
sol.par.g2_fun_arg = [0 1 1e-9 90];
sol.par.int2 = [0.1];
sol.par.g1_fun_arg = [0];
sol.par.tmax = 8e-10;
sol.par.tmesh_type = 'linear';
sol.par.side = 'right';
sol.par = refresh_device(sol.par);
sol_lightpulse = df(sol);

%%
delays = [10e-9 25e-9 50e-9 75e-9 100e-9 200e-9 400e-9 750e-9 1e-6];
J_pulse = cell(1,length(delays));
J_nopulse = cell(1,length(delays));

for i = 1:length(delays)
    display(num2str(i))
    delay = delays(i);
    sol_temp = sol_lightpulse;
    sol_temp.par.mobseti = 0;
    sol_temp.par.AbsTol = 1e-9;
    sol_temp.par.RelTol = 1e-6;
    sol_temp.par.g2_fun_type = 'constant';
    sol_temp.par.g2_fun_arg = [0];
    sol_temp.par.int2 = [0];
    sol_temp.par.tmax = delay;
    sol_temp.par = refresh_device(sol_temp.par);
    sol_delay = df(sol_temp);
    
    period = 2e-8;
    %As a percentage
    duty = 50;
    sol_delay.par.tmax = 0.95*(period);
    sol_delay.par.tpoints = 1e3;
    sol_nopulse = df(sol_delay);
    J_nopulse{i} = dfana.calcJ(sol_nopulse);

    sol_delay.par.PPP = 1;
    sol_delay.par.PPP_args = [0 1e5 period duty];
    sol_delay.par = refresh_device(sol_delay.par);
    
    sol_pulse = df(sol_delay);
    
    J_pulse{i} = dfana.calcJ(sol_pulse);
end

%%
figure('Name', 'PPP Current')

hold on
for i = 1:length(delays)
    plot(sol_pulse.t, (J_pulse{i}.tot(:,42)-J_nopulse{i}.tot(:,42))/J_lighton)
end
plot(sol_darkpulse.t, J_pulse_dark.tot(:,42)/J_lighton, 'Color', 'k')
hold off
xlabel('Time (s)')
ylabel('\DeltaJ/J')

%%
figure('Name', 'PPP Current')

hold on
for i = 1:length(delays)
    plot(sol_pulse.t, (J_pulse{i}.tot(:,42))/J_lighton)
end
plot(sol_darkpulse.t, J_pulse_dark.tot(:,42)/J_lighton, 'Color', 'k')
hold off
xlabel('Time (s)')
ylabel('\DeltaJ/J')

