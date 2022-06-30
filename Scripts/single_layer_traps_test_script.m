% Single layer test script
%
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
par = pc('Input_files/Kirchartz2011.csv');
%% No traps
par_notrap = par;
par_notrap.Nt_CB = 0;
par_notrap.Nt_VB = 0;
par_notrap.Cn_CB = 1e-100;
par_notrap.Cp_CB = 1e-100;
par_notrap.Cn_VB = 1e-100;
par_notrap.Cp_VB = 1e-100;
par_notrap = refresh_device(par_notrap);
% Load parameters

% Equilibrium oslutions
soleq_notrap = equilibrate(par_notrap, 1);

% Do JV scans
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV_notrap_dark = doCV(soleq_notrap.el, 0, 0, 1.0, 0, 1e-3, 1, 201);

sol_CV_notrap_1sun = doCV(soleq_notrap.el, 1, 0, 1.2, 0, 1e-3, 1, 201);

%% Traps
% Load parameters
par = pc('Input_files/Kirchartz2011.csv');
par_trap = par;
par_trap.Nt_CB = 1e19;
par_trap.Nt_VB = 1e19;
par_trap.Cn_CB = 1.1e-11;
par_trap.Cp_CB = 2.3e-10;
par_trap.Cn_VB = 5.2e-13;
par_trap.Cp_VB = 2.6e-10;
% par_trap.Cn_CB = 1e-10;
% par_trap.Cp_CB = 1e-10;
% par_trap.Cn_VB = 1e-10;
% par_trap.Cp_VB = 1e-10;
par_trap = refresh_device(par_trap);

% Equilibrium oslutions
soleq_trap = equilibrate(par_trap, 1);

% Do JV scans
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV_trap_dark = doCV(soleq_trap.el, 0, 0, 1.0, 0, 1e-3, 1, 201);

sol_CV_trap_1sun = doCV(soleq_trap.el, 1, 0, 1.2, 0, 1e-3, 1, 201);

%% Equilibrium np plots
dfplot.npx(soleq_notrap.el)
hold on
dfplot.npx(soleq_trap.el)
hold off

%% Equilibrium EL plots
dfplot.ELx(soleq_notrap.el)
hold on
dfplot.ELx(soleq_trap.el)
hold off

%% plot JV scan
dfplot.JtotVapp(sol_CV_notrap_dark, 0);
hold on
dfplot.JtotVapp(sol_CV_notrap_1sun, 0);
hold on
dfplot.JtotVapp(sol_CV_trap_dark, 0);
hold on
dfplot.JtotVapp(sol_CV_trap_1sun, 0);
hold off
legend('Nc0 = Nv0 = 0 cm-3, dark', 'Nc0 =  Nv0 = 0 cm-3, 1sun',...
    ['Nc0 = Nv0 =', num2str(par_trap.Nt_CB) ' cm-3, dark'],...
    ['Nc0 = Nv0 =', num2str(par_trap.Nt_CB) ' cm-3, 1 sun'])
ylim([-15e-3, 10e-3])

%% Fit dark curves
Vapp = dfana.calcVapp(sol_CV_notrap_dark);
J_notrap = dfana.calcJ(sol_CV_notrap_dark);
Jtot_notrap = J_notrap.tot(:,1);
fit_notrap = fit(Vapp(41:61)', log(Jtot_notrap(41:61)), 'poly1');
nid_notrap = (1/(par.kB*par.T))*(1/fit_notrap.p1);

Vapp = dfana.calcVapp(sol_CV_trap_dark);
J_trap = dfana.calcJ(sol_CV_trap_dark);
Jtot_trap = J_trap.tot(:,1);
fit_trap = fit(Vapp(16:36)', log(Jtot_trap(16:36)), 'poly1');
nid_trap = (1/(par.kB*par.T))*(1/fit_trap.p1);

%% Plot fits
figure(401)
plot(fit_notrap, 'k--', Vapp, log(Jtot_notrap), '-')
hold on
plot(fit_trap, 'k-.', Vapp, log(Jtot_trap), '-')
hold off
xlabel('Applied voltage (V)')
ylabel('ln(J)')
legend('Nc0 =  Nv0 = 0 cm-3, dark', 'fit, nid = 1.0', ...
    ['Nc0 = Nv0 =', num2str(par_trap.Nt_CB) ' cm-3, dark'], 'fit, nid = 1.4')
xlim([0,0.8])
ylim([-20,0])

%% Get trap densities
[n_trap_dark, ~] = dfana.calcn_trap(sol_CV_trap_dark, 'whole');
[p_trap_dark, ~] = dfana.calcp_trap(sol_CV_trap_dark, 'whole');
[n_trap_1sun, ~] = dfana.calcn_trap(sol_CV_trap_1sun, 'whole');
[p_trap_1sun, ~] = dfana.calcp_trap(sol_CV_trap_1sun, 'whole');

%% Make movies
% ntrap vs position vs energy
n_trap_Et_movie(sol_CV_trap_dark, 'n_trap_1_layer', [-1.2, 0], [1e10, 1e17], 1);

%%
makemovie(sol_CV_trap_dark, @dfplot.ELnpx_trap, 0, [1e12, 1e18], 'ELnpx_1sun', 1, 0, n_trap_dark, p_trap_dark);
makemovie(sol_CV_trap_1sun, @dfplot.ELnp_trap, 0, [1e12, 1e18], 'ELnpx_1sun', 1, 0, n_trap_1sun, p_trap_1sun);
