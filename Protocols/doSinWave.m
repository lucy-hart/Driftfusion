function [sol_pulse] = doSinWave(sol_ini, pulse_int, tmax, tpoints, frequency, mobseti)
% Uses square wave light generator for light source 2 and peforms a single
% pulse
%
%% Input arguments
% SOL_INI = initial conditions
% PULSE_INT
% MOBSETI = Ion mobility switch
% RS = Series resistance - recommended to use Rs = 1e6 for approx open
% circuit
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
disp('Starting')
par = sol_ini.par;


% Setup time mesh
par.tmesh_type = 1;
par.tmax = tmax;
par.t0 = 0;

par.mobseti = mobseti;

%% Setup square wave function generator
par.g2_fun_type = 'sin';
par.tpoints = tpoints;
par.g2_fun_arg(1) = 0.5*pulse_int;          % Offset
par.g2_fun_arg(2) = 0.5*pulse_int;          % Amplitude of sin wave
par.g2_fun_arg(3) = frequency;              % Frequency
par.g2_fun_arg(4) = 0;                      % Phase Shift
sol_pulse = df(sol_ini, par);

disp('Complete')

end