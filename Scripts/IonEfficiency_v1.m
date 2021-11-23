%Program to do a parameter sweep of E_bi and recombination velocity for a single
%layer MAPI cell and see how device parameters are affected by the presence
%of ions.

%% Read in basic device
par=pc('Input_files/1_layer_test.csv');

%% Define parameter space
% Rows are recombination velocities, columns are E_bi values
params = zeros(7,12);
sigma_values = logspace(1,7,7);
E_BI_values = linspace(0,1.2,12);
for i=1:7
    for j=1:12
        params(i,j) = [sigma_values()]

%% Do (many) JV sweeps

%% Plot results