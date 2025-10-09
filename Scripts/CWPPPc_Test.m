par = pc('Input_files/SAM_MAFACsPbIBr_C60.csv');
par.vsr_mode = 0;
par.Rs = 0;
par.z_t = -1;
par = refresh_device(par);
%%
eqm = equilibrate(par);
eqm.par.RelTol = 1e-5;
eqm.par.AbsTol = 1e-8;

%%
powers_pump = linspace(1,10,5);
powers_push = linspace(10,70,7);
J_lighton = zeros(length(powers_pump),1);
DeltaJ_IR = zeros(length(powers_push), length(powers_pump));

%%
for j = 1:length(powers_pump)
    sol = eqm.ion;
    sol.par.laser_lambda2 = 750;
    sol.par.g2_fun_type = 'constant';
    sol.par.g2_fun_arg = [1];
    sol.par.int2 = [powers_pump(j)];
    sol.par.g1_fun_arg = [0];
    sol.par.tmax = 1;
    sol.par.tmesh_type = 'linear';
    sol.par.side = 'right';
    sol.par = refresh_device(sol.par);
    sol_light = df(sol);
    

    J = dfana.calcJ(sol_light);
    t_start_mean = 0.5;
    t_end_mean = 0.9;
    [~,arg_start_mean] = min(abs(sol_light.t - t_start_mean));
    [~,arg_end_mean] = min(abs(sol_light.t - t_end_mean));
    
    J_lighton(j) = mean(J.tot(arg_start_mean:arg_end_mean,1));
    

    for i = 1:length(powers_push)
        sol = eqm.ion;
        sol.par.g2_fun_arg = [0];
        sol.par.int2 = [0];
        sol.par.g1_fun_arg = [0];
        sol.par.tmesh_type = 'linear';
        sol.par.side = 'right';
        
        period = 40;
        %As a percentage
        duty = 50;
        power = powers_push(i);
        tstop = 0.99*period;
        tpoints = 500;
        sol.par.tmax = tstop;
        sol.par.tpoints = tpoints;
        sol.par.PPP = 1;
        sol.par.PPP_args = [0 power period duty];
        sol.par = refresh_device(sol.par);
        
        sol_darkpulse = df(sol);
        J_pulse_dark = dfana.calcJ(sol_darkpulse).tot;
            
        sol = sol_light;
        
        sol.par.tmax = tstop;
        sol.par.tpoints = tpoints;
        sol.par.PPP = 1;
        sol.par.PPP_args = [0 power period duty];
        sol.par = refresh_device(sol.par);
        
        sol_lightpulse = df(sol);
        J_pulse_light = dfana.calcJ(sol_lightpulse).tot;

        Jnet = J_pulse_light(:,42) - J_pulse_dark(:,42);

        t_start_mean = 0.2;
        t_end_mean = 0.4;
        [~,arg_start_mean] = min(abs(sol_lightpulse.t - t_start_mean*period));
        [~,arg_end_mean] = min(abs(sol_lightpulse.t - t_end_mean*period));
        DeltaJ_IR(i,j) = mean(Jnet(arg_start_mean:arg_end_mean)) - J_lighton(j);
    end
end

%%
figure('Name', 'PPP Current')
colours = ["red", "blue", "green", "black", "yellow"];
hold on
for j = 1:length(powers_pump)
    for i = 1:length(powers_push)
        plot(powers_push, DeltaJ_IR(:,j)/J_lighton(j), 'Marker', 'o', 'Color', colours(j))
    end
end
hold off
xlabel('IR Push Strength (Arb.)')
ylabel('\DeltaJ/J')
%ylim([4e-5 9e-5])