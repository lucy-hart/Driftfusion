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
par = pc('Input_files/pn_heterojunction.csv');
%% No traps
par_notrap = par;
par_notrap.Nt_CB = [0, 0, 0];
par_notrap.Nt_VB = [0, 0, 0];
par_notrap.taun_CB = [1e6, 1e20, 1e6];
par_notrap.taup_CB = [1e6, 1e20, 1e6];
par_notrap.taun_VB = [1e6, 1e20, 1e6];
par_notrap.taup_VB = [1e6, 1e20, 1e6];
par_notrap = refresh_device(par_notrap);
% Load parameters

% Equilibrium solution
soleq_notrap = equilibrate(par_notrap, 1);

% Do JV scans
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV_notrap = doCV(soleq_notrap.el, 1, 0, 0.8, 0, 100e-3, 1, 201);

%% Traps
% Load parameters
par_trap = par;
par_trap.Nt_CB = [0, 0, 1e16];
par_trap.Nt_VB = [1e14, 0, 0];
par_trap.taun_CB = [1e20, 1e20, 1e-6];
par_trap.taup_CB = [1e20, 1e20, 1e-6];
par_trap.taun_VB = [1e-6, 1e20, 1e20];
par_trap.taup_VB = [1e-6, 1e20, 1e20];
par_trap = refresh_device(par_trap);

% Equilibrium solution
soleq_trap = equilibrate(par_trap, 1);

% Do JV scans
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV_trap = doCV(soleq_trap.el, 1, 0, 0.8, 0, 100e-3, 1, 201);

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
dfplot.JtotVapp(sol_CV_notrap, 0);
hold on
dfplot.JtotVapp(sol_CV_trap, 0);
hold off
legend('no traps', 'traps')
ylim([-20e-3, 10e-3])
%set(gca,'YScale','log')