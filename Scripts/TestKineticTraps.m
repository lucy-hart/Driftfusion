par = pc('Input_files/SAM_MAFACsPbIBr_C60.csv');
par.vsr_mode = 0;
par.Rs = 1e6;
par = refresh_device(par);
eqm = equilibrate(par,0,1);
doJV = 1;
%%
if doJV == 1
    JVsol = doCV(eqm.ion, 1, -0.2, 1.2, -0.2, 1e-3, 1, 281);
    %Plot_Current_Contributions(JVsol)
end

%%
t_stop = 1e-3;
%As a percentage
duty = 0.00005;
power = 0.1;

eqm.ionstatic.par.laser_lambda2 = 605;
eqm.ionstatic.par.g2_fun_type = 'square';
eqm.ionstatic.par.g2_fun_arg = [0 power t_stop duty];
eqm.ionstatic.par.g1_fun_arg = [0];
eqm.ionstatic.par.tmax = 0.95*t_stop;
eqm.ionstatic.par.tmesh_type = 'linear';
eqm.ionstatic.par.tpoints = 5e3;
eqm.ionstatic.par.side = 'right';
eqm.ionstatic.par = refresh_device(eqm.ionstatic.par);

pulse_sol = df(eqm.ionstatic);

t_eqm = pulse_sol.t;
PL_eqm = dfana.calcPLt(pulse_sol);

eqm.ion.par.laser_lambda2 = 605;
eqm.ion.par.g2_fun_type = 'square';
eqm.ion.par.g2_fun_arg = [0 power t_stop duty];
eqm.ion.par.g1_fun_arg = [0];
eqm.ion.par.tmax = 0.95*t_stop;
eqm.ion.par.tmesh_type = 'linear';
eqm.ion.par.tpoints = 5e3;
eqm.ion.par.kineticset = 1;
eqm.ion.par.side = 'right';
eqm.ion.par = refresh_device(eqm.ion.par);

pulse_sol_kinetic = df(eqm.ion);

t_kinetic = pulse_sol_kinetic.t;
PL_kinetic = dfana.calcPLt(pulse_sol_kinetic);
[~,argstart_k] = min(abs(t_kinetic-t_stop*duty/100));
[~,argstart_e] = min(abs(t_eqm-t_stop*duty/100));
%%
figure('Name', 'ComparePL')

hold on
plot(t_eqm(argstart_e:end), PL_eqm(argstart_e:end)/max(PL_eqm(argstart_e:end)), 'Color', 'Red')
plot(t_kinetic(argstart_k:end), PL_kinetic(argstart_k:end)/max(PL_kinetic(argstart_k:end)), 'Color', 'Blue')
hold off

xlabel('Time (s)')
ylabel('PL (norm.)')
ylim([5e-4 1.1])
xlim([0 t_stop])
set(gca, 'YScale', 'log')
%set(gca, 'XScale', 'log')

%%
figure('Name', 'ComparePL')

hold on
plot(t_eqm, PL_eqm, 'Color', 'Red')
plot(t_kinetic, PL_kinetic, 'Color', 'Blue')
hold off

xlabel('Time (s)')
ylabel('PL')
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')


%%
figure('Name', 'ComparePL')

hold on
plot(t_eqm, PL_eqm./PL_kinetic, 'Color', 'Blue')
%plot(t_kinetic, PL_kinetic, 'Color', 'Blue')
hold off

xlabel('Time (s)')
ylabel('PL')
ylim([0 1.1])
%set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')

