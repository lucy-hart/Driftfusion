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
%% Start code
% Load parameters
par = pc('Input_files/1_layer_test.csv');

%% Equilibrium oslutions
soleq = equilibrate(par);
dfplot.npx(soleq.el)

%% Do JV scans
% sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)
sol_CV = doCV(soleq.el, 1, 0, 1, 0, 100e-3, 1, 201);

%% plot JV scan
dfplot.JVapp(sol_CV, 0);

%set(gca,'YScale','log')