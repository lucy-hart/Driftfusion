%Program to do a parameter sweep of Ebi and recombination velocity for a single
%layer MAPI cell and see how device parameters are affected by the presence
%of ions.


%% Define parameter space
% Rows are recombination velocities, columns are Ebi values
params = cell(7,12);
sigma_values = logspace(1,7,7);
Ebi_values = linspace(0,0.6,12);
for i=1:7
    for j=1:12
        params{i,j} = [sigma_values(i), Ebi_values(j)];
    end
end

%% Read in basic device
par=pc('Input_files/1_layer_test.csv');

%% Do (many) JV sweeps
ion_results = cell(7,12);
el_results = cell(7,12);
for i=1:7
    for j=1:12
        par.Phi_left = par.EF0-params{i,j}(2);
        par.Phi_right = par.EF0+params{i,j}(2);
        par.sn_l = params{i,j}(2);
        par.sp_r = params{i,j}(2);
        par = refresh_device(par);
        soleq = equilibrate(par);
        solCV_ion = doCV(soleq.ion, 1, -0.1, 2*params{i,j}(2) + 0.4, -0.1, 10e-3, 1, 241);
        solCV_el = doCV(soleq.el, 1, -0.1, 2*params{i,j}(2) + 0.4, -0.1, 10e-3, 1, 241);
        ion_results{i,j} = CVstats(solCV_ion);
        el_results{i,j} = CVstats(solCV_el);
    end
end

%% Plot results