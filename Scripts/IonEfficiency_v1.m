%Program to do a parameter sweep of Vbi and recombination velocity for a single
%layer MAPI cell and see how device parameters are affected by the presence
%of ions.
tic
%% Define parameter space
% Rows are recombination velocities, columns are Ebi values
params = cell(7,13);
sigma_values = logspace(1,7,7);
Ebi_values = zeros(1,14);
Ebi_values(2:end) = linspace(0,0.6,13);
Ebi_values(2) = 0.01;

for i=1:7
    for j=1:14
            params{i,j} = [sigma_values(i), Ebi_values(j)];
    end
end

ion_results = cell(7,14);
el_results = cell(7,14);
%% Do (many) JV sweeps
par=pc('Input_files/1_layer_test.csv');
for i=1:7
    for j=1:14
        if j==1 && i ==7
            par.Phi_left = par.EF0-0.01;
            par.Phi_right = par.EF0+0.01;
        else
            par.Phi_left = par.EF0-params{i,j}(2);
            par.Phi_right = par.EF0+params{i,j}(2);  
        end
        par.sn_l = params{i,j}(1);
        par.sp_r = params{i,j}(1);
        par = refresh_device(par);
        soleq = equilibrate(par);
        if params{i,j}(2) <= 0.35
            solCV_ion = doCV(soleq.ion, 1, -0.25, 1.1, -0.25, 1e-4, 1, 241);
            solCV_el = doCV(soleq.el, 1, -0.25, 1.1, -0.25, 1e-4, 1, 241);
        elseif params{i,j}(2) > 0.35
            solCV_ion = doCV(soleq.ion, 1, -0.25, 1.25, -0.25, 1e-4, 1, 241);
            solCV_el = doCV(soleq.el, 1, -0.25, 1.25, -0.25, 1e-4, 1, 241);
        end
        ion_results{i,j} = CVstats(solCV_ion);
        el_results{i,j} = CVstats(solCV_el);
    end
end

toc
%% Calculate figure of merit
PCE_ratio = zeros(7,14);
for i=1:7
    for j=1:14
        try
            PCE_ratio(i,j) = el_results{i,j}.efficiency_f/ion_results{i,j}.efficiency_f;
        catch
            warning('CVstats unsucessful. Assigning a PCE value of 0.');
            PCE_ratio(i,j) = 0;
        end
    end
end

%% Plot results 
PCE_ratio(6:7,1) = 0;
figure(1)
contourf(2*Ebi_values, log10(sigma_values), log10(PCE_ratio), 12, 'LineWidth', 0.1)
xlabel('V_{BI} (V)')
ylabel('log_{10}(s_{surf} /cms^{-1})')
c = colorbar;
c.Label.String = 'log_{10}(PCE_{el} / PCE_{ion})';

%% Voc vs V_bi plot

V_ratio = zeros(7,13);
for i = 1:7
    for j = 1:13
        V_ratio(i,j) = ion_results{i,j+1}.Voc_f/(2*Ebi_values(j+1));
    end
end
figure(2)
contourf(2*Ebi_values(2:end), log10(sigma_values), log10(V_ratio), 12, 'LineWidth', 1)
xlabel('V_{BI} (V)')
ylabel('log_{10}(s_{surf} /cms^{-1})')
c = colorbar;
c.Label.String = 'log_{10}(V_{OC,ion} / V_{BI})';

