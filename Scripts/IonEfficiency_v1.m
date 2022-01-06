%Program to do a parameter sweep of Vbi and recombination velocity for a single
%layer MAPI cell and see how device parameters are affected by the presence
%of ions.
tic
%% Define parameter space
% Rows are recombination velocities, columns are Ebi values
sigma_values = logspace(1, 4, 5);
Ebi_values = linspace(0, 0.5, 6);
Nsigma = length(sigma_values); 
NEbi = length(Ebi_values);
params = cell(Nsigma, NEbi);

for i=1:Nsigma
    for j=1:NEbi
            params{i,j} = [sigma_values(i), Ebi_values(j)];
    end
end

error_log = zeros(Nsigma, NEbi);
ion_results = cell(Nsigma, NEbi);
el_results = cell(Nsigma, NEbi);
%% Do (many) JV sweeps
par=pc('Input_files/1_layer_test.csv');
for i=1:Nsigma
    for j=1:NEbi
        disp(["minority carrier SRV = ", num2str(sigma_values(i)), "cm s-1"])
        disp(["Built in potential, Vbi = ", num2str(Ebi_values(j)), " V"])
%         if j== 1 && i == Nsigma
%             par.Phi_left = par.EF0 - 0.01;
%             par.Phi_right = par.EF0 + 0.01;
%         else
            par.Phi_left = par.EF0 - params{i,j}(2);
            par.Phi_right = par.EF0 + params{i,j}(2);  
%         end
        par.sn_l = params{i,j}(1);
        par.sp_r = params{i,j}(1);
        par = refresh_device(par);
        
        soleq{i, j} = equilibrate(par);
        try
            if params{i,j}(2) <= 0.35
                solCV_ion{i, j} = doCV(soleq{i, j}.ion, 1, -0.1, 1, -0.1, 1e-4, 1, 241);
                solCV_el{i, j} = doCV(soleq{i, j}.el, 1, -0.1, 1, -0.1, 1e-4, 1, 241);
            elseif params{i,j}(2) > 0.35
                solCV_ion{i, j} = doCV(soleq{i, j}.ion, 1, -0.1, 1.2, -0.1, 1e-4, 1, 241);
                solCV_el{i, j} = doCV(soleq{i, j}.el, 1, -0.1, 1.2, -0.1, 1e-4, 1, 241);
            end
            error_log(i,j) = 0;
            ion_results{i,j} = CVstats(solCV_ion{i, j});
            el_results{i,j} = CVstats(solCV_el{i, j});
        catch
            warning("JV solution failed, try reducing Vmax")
            error_log(i,j) = 1;
        end
    end
end

toc
%% Calculate figure of merit
PCE_ratio = zeros(Nsigma, NEbi);
for i=1:Nsigma
    for j=1:NEbi
        try
            PCE_ratio(i,j) = el_results{i,j}.efficiency_f/ion_results{i,j}.efficiency_f;
            PCE_el(i,j) = el_results{i,j}.efficiency_f;
            PCE_ion(i,j) = ion_results{i,j}.efficiency_f;
        catch
            warning('CVstats unsucessful. Assigning a PCE value of 0.');
            PCE_el(i,j) = NaN;
            PCE_ion(i,j) = NaN;
            PCE_ratio(i,j) = NaN;
        end
    end
end

%% Plot results 
figure(1)
contourf(2*Ebi_values, log10(sigma_values), PCE_ratio, 12, 'LineWidth', 0.1)
xlabel('V_{BI} (V)')
ylabel('log_{10}(s_{surf} /cms^{-1})')
c = colorbar;
c.Label.String = 'PCE_{el} / PCE_{ion}';

%% Plot results 
figure(12)
contourf(2*Ebi_values, log10(sigma_values), PCE_ion, 12, 'LineWidth', 0.1)
xlabel('V_{BI} (V)')
ylabel('log_{10}(s_{surf} /cms^{-1})')
c = colorbar;
c.Label.String = 'PCE_{ion}';

%% Plot results 
figure(13)
contourf(2*Ebi_values, log10(sigma_values), PCE_el, 12, 'LineWidth', 0.1)
xlabel('V_{BI} (V)')
ylabel('log_{10}(s_{surf} /cms^{-1})')
c = colorbar;
c.Label.String = 'PCE_{el}';

%% Voc vs V_bi plot

V_ratio = zeros(Nsigma,NEbi);
for i = 1:Nsigma
    for j = 1:NEbi
        V_ratio(i,j) = ion_results{i,j}.Voc_f/(2*Ebi_values(j));
    end
end
figure(2)
contourf(2*Ebi_values(1:end), log10(sigma_values), log10(V_ratio), 12, 'LineWidth', 1)
xlabel('V_{BI} (V)')
ylabel('log_{10}(s_{surf} /cms^{-1})')
c = colorbar;
c.Label.String = 'log_{10}(V_{OC,ion} / V_{BI})';

