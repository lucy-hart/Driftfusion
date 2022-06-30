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
par_PSC = pc('Input_files/spiro_mapi_tio2.csv');
%% No traps
par_notrap_PSC = par_PSC;
par_notrap_PSC.E_UCB = [0.07, 0.07, 0.07, 0.07, 0.07];
par_notrap_PSC.E_UVB = [0.07, 0.07, 0.07, 0.07, 0.07];
par_notrap_PSC.Nt_CB = [0,0,0,0,0];
par_notrap_PSC.Nt_VB = [0,0,0,0,0];
par_notrap_PSC.Cn_CB = [1e-100,1e-100,1e-100,1e-100,1e-100];
par_notrap_PSC.Cp_CB = [1e-100,1e-100,1e-100,1e-100,1e-100];
par_notrap_PSC.Cn_VB = [1e-100,1e-100,1e-100,1e-100,1e-100];
par_notrap_PSC.Cp_VB = [1e-100,1e-100,1e-100,1e-100,1e-100];
par_notrap_PSC = refresh_device(par_notrap_PSC);
% Load parameters

% Equilibrium oslutions
soleq_notrap_PSC = equilibrate(par_notrap_PSC);

% Do JV scans
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV_notrap_PSC_dark = doCV(soleq_notrap_PSC.ion, 0, 0, 1.2, 0, 100e-3, 1, 201);

sol_CV_notrap_PSC = doCV(soleq_notrap_PSC.ion, 1, 0, 1.2, 0, 100e-3, 1, 201);

%% Traps
% Load parameters
par_trap_PSC = par_PSC;
par_trap_PSC.E_UCB = [0.07, 0.07, 0.07, 0.07, 0.07];
par_trap_PSC.E_UVB = [0.07, 0.07, 0.07, 0.07, 0.07];
par_trap_PSC.Nt_CB = [1e10,0,1e14,0,1e20];
par_trap_PSC.Nt_VB = [1e20,0,1e16,0,1e10];
par_trap_PSC.Cn_CB = [1e-12,1e-100,1e-100,1e-100,1e-12];
par_trap_PSC.Cp_CB = [1e-10,1e-100,1e-100,1e-100,1e-11];
par_trap_PSC.Cn_VB = [1e-11,1e-100,1e-100,1e-100,1e-10];
par_trap_PSC.Cp_VB = [1e-12,1e-100,1e-100,1e-100,1e-12];
% par_trap_PSC.Cn_CB = 1e-10;
% par_trap_PSC.Cp_CB = 1e-10;
% par_trap_PSC.Cn_VB = 1e-10;
% par_trap_PSC.Cp_VB = 1e-10;
par_trap_PSC = refresh_device(par_trap_PSC);

% Equilibrium oslutions
soleq_trap_PSC = equilibrate(par_trap_PSC);

% Do JV scans
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV_trap_PSC_dark = doCV(soleq_trap_PSC.ion, 0, 0, 1.2, 0, 100e-3, 1, 201);

sol_CV_trap_PSC = doCV(soleq_trap_PSC.ion, 1, 0, 1.2, 0, 100e-3, 1, 201);

%% Equilibrium np plots
dfplot.npx(soleq_notrap_PSC.ion)
hold on
dfplot.npx(soleq_trap_PSC.ion)
hold off

%% Equilibrium EL plots
dfplot.ELx(soleq_notrap_PSC.ion)
hold on
dfplot.ELx(soleq_trap_PSC.ion)
hold off

%% plot JV scan
dfplot.JtotVapp(sol_CV_notrap_PSC_dark, 0);
hold on
dfplot.JtotVapp(sol_CV_notrap_PSC, 0);
hold on
dfplot.JtotVapp(sol_CV_trap_PSC_dark, 0);
hold on
dfplot.JtotVapp(sol_CV_trap_PSC, 0);
legend('no traps, dark', 'no traps, 1 sun', 'traps, dark', 'traps, 1 sun')
hold off
ylim([-30e-3, 10e-3])

%% Get trap densities
[n_trap_PSC, ~] = dfana.calcn_trap(sol_CV_trap_PSC, 'whole');
[p_trap_PSC, ~] = dfana.calcp_trap(sol_CV_trap_PSC, 'whole');

%% Make movies
ELnpacx_1sun = makemovie_trap(sol_CV_trap_PSC, @dfplot.ELxnpxacx_trap, 0, 0, 'ELxnpxacx_1sun', 1, 0, n_trap_PSC, p_trap_PSC);

%% ntrap vs position vs energy
n_trap_Et_movie(sol_CV_trap_PSC_dark, 'n_trap_PSC', [-6, -3], [1, 1e17], 1);