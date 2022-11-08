%Use this file for energy related things (TL Fermi levels, energy offsets
%and Vbi)

tic
%% Define parameter space
%Rows are HOMO offsets
%Columns are LUMO offsets

n_values = 11;
Delta_HOMO = linspace(-0.2, 0.3, n_values);
Delta_LUMO = linspace(-0.3, 0.2, n_values);

params = cell(n_values, n_values);

for i=1:n_values
    for j=1:n_values
            params{i,j} = [Delta_HOMO(i), Delta_LUMO(j)];
    end
end

%%
error_log_el = zeros(n_values, n_values);
error_log_ion = zeros(n_values, n_values);
soleq = cell(n_values, n_values);
solCV_ion = cell(n_values, n_values);
solCV_el = cell(n_values, n_values);
ion_results = cell(n_values, n_values);
el_results = cell(n_values, n_values);

%% Do (many) JV sweeps
par=pc('Input_files/EnergyOffsetSweepParameters.csv');

%Reset the electrode work functions in each loop to be safe as they are
%changed for the cases where E_LUMO (E_HOMO) is far below (above) the CB
%(VB)
for i= 1:n_values
    for j= 1:n_values
        disp(["DeltaE_HOMO = ", num2str(Delta_HOMO(i)), " eV"])
        disp(["DeltaE_LUMO = ", num2str(Delta_LUMO(j)), " eV"])
        %HTL Energetics
        par.Phi_left = -5.3;
        par.Phi_IP(1) = par.Phi_IP(3) + params{i,j}(1);
        par.Phi_EA(1) = par.Phi_IP(1) + 2.5;
        par.EF0(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
        par.Et(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
        if par.Phi_left < par.Phi_IP(1)
            par.Phi_left = par.Phi_IP(1);
        end
        %ETL Energetics
        par.Phi_right = -4.1;
        par.Phi_EA(5) = par.Phi_EA(3) + params{i,j}(2);
        par.Phi_IP(5) = par.Phi_EA(5) - 2.5;
        par.EF0(5) = (par.Phi_IP(5)+par.Phi_EA(5))/2;
        par.Et(5) = (par.Phi_IP(5)+par.Phi_EA(5))/2;
        if par.Phi_right > par.Phi_EA(5)
            par.Phi_right = par.Phi_EA(5);
        end

        par = refresh_device(par);

        soleq{i,j} = equilibrate(par);
        try            
            if params{i,j}(1) > 0.15 || params{i,j}(2) < - 0.15
                solCV_el{i, j} = doCV(soleq{i, j}.el, 1, -0.2, 1.0, -0.2, 1e-4, 1, 241);
            elseif (0 < params{i,j}(1)) && (params{i,j}(1)< 0.15) || (params{i,j}(2) > - 0.15) && (params{i,j}(2) < 0)
                solCV_el{i, j} = doCV(soleq{i, j}.el, 1, -0.2, 1.15, -0.2, 1e-4, 1, 271);
            else
                solCV_el{i, j} = doCV(soleq{i, j}.el, 1, -0.2, 1.2, -0.2, 1e-4, 1, 281);
            end
            error_log_el(i,j) = 0;
            el_results{i,j} = CVstats(solCV_el{i, j});
        catch
            warning("Electronic-only JV solution failed, try reducing Vmax")
            error_log_el(i,j) = 1;
        end

        try            
            if params{i,j}(1) > 0.15 || params{i,j}(2) < - 0.15
                solCV_ion{i, j} = doCV(soleq{i, j}.ion, 1, -0.2, 1.1, -0.2, 1e-4, 1, 271);
            else
                solCV_ion{i, j} = doCV(soleq{i, j}.ion, 1, -0.2, 1.2, -0.2, 1e-4, 1, 281);
            end
            error_log_ion(i,j) = 0;
            ion_results{i,j} = CVstats(solCV_ion{i, j});
        catch
            warning("JV solution failed, try reducing Vmax")
            error_log_ion(i,j) = 1;
        end
    end
end

toc

%% Plot results 
Voc_array_ion = zeros(n_values, n_values);
Voc_array_el = zeros(n_values, n_values);
for i = 1:n_values
    for j = 1:n_values
        try
            Voc_array_ion(i,j) = ion_results{i,j}.Voc_f;
        catch
            warning('No Voc value')
            Voc_array_ion(i,j) = 0;
        end

        try
            Voc_array_el(i,j) = el_results{i,j}.Voc_f;
        catch
            warning('No Voc value')
            Voc_array_el(i,j) = 0;
        end
    end
end 

%%
figure(1)
contourf(Delta_HOMO, Delta_LUMO, Voc_array_el, [0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20], 'LineWidth', 0.1)
xlabel('\DeltaE_{HOMO} (eV)')
ylabel('\DeltaE_{LUMO} (eV)')
ylim([-0.3, 0.2])
c = colorbar;
c.Label.String = 'V_{OC}';

%%
figure(6)
LevelList = linspace(-23,0,24);

contourf(voltage_matrix_plot(1:end,2:120)', ...
    (2*y_values(1:end)'*ones(1,119))', ...
    1000*result_matrix(1:end,2:120)', ...
    'LevelList', LevelList, ...
    'LineStyle', 'none')
hold on
contour(voltage_matrix_plot(1:end,2:120)', ...
       (2*y_values(1:end)'*ones(1,119))', ...
       1000*result_matrix(1:end,2:120)', [0 0], ...
       'LineWidth', 2, 'LineStyle', '-', 'LineColor', 'k', 'ShowText', 'off') 
plot(2*y_values(1:14), 2*y_values(1:14), 'k--')
hold off

cmap = parula;
%Reverse colourmap
%c = flipud(c);
colormap(cmap)
ylabel('V_{BI} (V)')
xlabel('Applied Voltage (V)')
ylim([0,1.2])
xlim([-0.19,1.19])
c = colorbar;
c.Label.String = 'Current Density (mAcm^{-3})';

%% Save results and solutions

%filename = 'tld_symmetric_Vbi_vs_Ncat_Piers_vsr.mat';
%save(filename, 'el_results', 'ion_results', 'solCV_el', 'solCV_ion')
