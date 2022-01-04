%Program to do a parameter sweep of Ebi and recombination velocity for a single
%layer MAPI cell and see how device parameters are affected by the presence
%of ions.

%% Define parameter space
% Rows are recombination velocities, columns are Ebi values
params = cell(7,13);
sigma_values = logspace(1,7,7);
Ebi_values = linspace(0,0.6,13);
Ebi_values(1) = 0.01;
for i=1:7
    for j=1:13
        params{i,j} = [sigma_values(i), Ebi_values(j)];
    end
end

ion_results = cell(7,13);
el_results = cell(7,13);
%% Do (many) JV sweeps
par=pc('Input_files/1_layer_test.csv');
for i=1:7
    for j=1:13
        par.Phi_left = par.EF0-params{i,j}(2);
        par.Phi_right = par.EF0+params{i,j}(2);  
        par.sn_l = params{i,j}(1);
        par.sp_r = params{i,j}(1);
        par = refresh_device(par);
        soleq = equilibrate(par);
        solCV_ion = doCV(soleq.ion, 1, -0.2, 1.1, -0.2, 1e-4, 1, 241);
        solCV_el = doCV(soleq.el, 1, -0.2, 1.1, -0.2, 1e-4, 1, 241);
        ion_results{i,j} = CVstats(solCV_ion);
        el_results{i,j} = CVstats(solCV_el);
    end
end

%% Plot results
PCE_ratio = zeros(7,12);
for i=1:7
    for j=1:13
        try
            PCE_ratio(i,j) = el_results{i,j}.efficiency_f/ion_results{i,j}.efficiency_f;
        catch
            warning('CVstats unsucessful. Assigning a PCE value of 0.');
            PCE_ratio(i,j) = 0;
        end
    end
end
figure(1)
contourf(2*Ebi_values, log10(sigma_values), log10(PCE_ratio), 10, 'LineWidth', 1)
xlabel('V_{BI} (V)')
ylabel('log_{10}(s_{surf} /cms^{-1})')
c = colorbar;
c.Label.String = 'log_{10}(PCE_{el} /PCE_{ion})';
