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
par=pc('Input_files/EnergyOffsetSweepParameters_v2.csv');

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

        Voc_max = 1.2;
        num_points = 281;
        while Voc_max >= 1.05
            try            
                solCV_el{i, j} = doCV(soleq{i, j}.el, 1, -0.2, Voc_max, -0.2, 1e-4, 1, num_points);           
                error_log_el(i,j) = 0;
                el_results{i,j} = CVstats(solCV_el{i, j});
                Voc_max = 0;                
            catch
                if Voc_max > 1.05
                    warning("Electronic-only JV solution failed, reducing Vmax by 0.05 V")
                    Voc_max = Voc_max - 0.03;
                    num_points = num_points - 6;
                elseif Voc_max == 1.05
                    warning("Electronic-only JV solution failed.")
                    error_log_el(i,j) = 1;
                    el_results{i,j} = 0;
                end
            end
        end

        Voc_max = 1.2;
        num_points = 281; 
        while Voc_max >= 1.05
            try
                solCV_ion{i, j} = doCV(soleq{i, j}.ion, 1, -0.2, Voc_max, -0.2, 1e-4, 1, num_points);
                error_log_ion(i,j) = 0;
                ion_results{i,j} = CVstats(solCV_ion{i, j});
                Voc_max = 0;
            catch
                if Voc_max > 1.05
                    warning("Ionic JV solution failed, reducing Vmax by 0.05 V")
                    Voc_max = Voc_max - 0.03;
                    num_points = num_points - 6;
                elseif Voc_max == 1.05
                    warning("Ionic JV solution failed.")
                    error_log_ion(i,j) = 1;
                    ion_results{i,j} = 0;
                end
            end
        end
    end
end

toc

%% Plot results 
Stats_array_ion = zeros(n_values, n_values, 4);
Stats_array_el = zeros(n_values, n_values, 4);
for i = 1:n_values
    for j = 1:n_values
        try
            Stats_array_ion(i,j,1) = 1e3*ion_results{i,j}.Jsc_f;
            Stats_array_ion(i,j,2) = ion_results{i,j}.Voc_f;
            Stats_array_ion(i,j,3) = ion_results{i,j}.FF_f;
            Stats_array_ion(i,j,4) = ion_results{i,j}.efficiency_f;
        catch
            warning('No Stats (ion)')
            Stats_array_ion(i,j,:) = 0;
        end

        try
            Stats_array_el(i,j,1) = 1e3*el_results{i,j}.Jsc_f;
            Stats_array_el(i,j,2) = el_results{i,j}.Voc_f;
            Stats_array_el(i,j,3) = el_results{i,j}.FF_f;
            Stats_array_el(i,j,4) = el_results{i,j}.efficiency_f;
        catch
            warning('No Stats (el)')
            Stats_array_el(i,j,:) = 0;
        end
    end
end 

%%
figure('Name', 'JV Params vs Energy Offsets')
ion = 0;
num = 1;
labels = ["J_{SC} (mA cm^{-2})", "V_{OC} (V)", "FF", "PCE (%)"];
lims = [[-23 -19]; [0.77 1.17]; [0.1 0.85]; [1 21]];
if ion == 1 
    data = Stats_array_ion;
elseif ion == 0
    data = Stats_array_el;
end
%Coutourf has coumn and row inidices as x and y respectively
contourf(Delta_LUMO, Delta_HOMO, data(:,:,num), 'LineWidth', 0.1)
xlabel('\DeltaE_{LUMO} (eV)')
ylabel('\DeltaE_{HOMO} (eV)')
ylim([-0.2, 0.3])
c = colorbar;
c.Label.String = labels(num);
clim(lims(num,:))

%%
figure('Name', 'Normalised PCE vs Energy Offsets')
num = 4;
labels = ["J_{SC} (mA cm^{-2})", "V_{OC} (V)", "FF", "$\mathrm{\frac{PCE_{ion}}{PCE_{no ion}}}$ - 1"];
if ion == 1 
    data = Stats_array_ion;
elseif ion == 0
    data = Stats_array_el;
end
%Coutourf has coumn and row inidices as x and y respectively
contourf(Delta_LUMO, Delta_HOMO, (Stats_array_ion(:,:,num) - Stats_array_el(:,:,num))./(Stats_array_el(:,:,num)), ...
    'SHowText', true, 'LineWidth', 0.1)
xlabel('\DeltaE_{LUMO} (eV)')
ylabel('\DeltaE_{HOMO} (eV)')
ylim([-0.2, 0.3])
c = colorbar;
c.Label.Interpreter = 'latex';
c.Label.String = labels(num);
c.Label.FontSize = 20;
clim([-0.85 0.05])

%%
doFigure = 0; 
if doFigure == 1
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
end



%% Save results and solutions
save_file = 0;
if save_file == 1
    filename = 'tld_symmetric_DeltaEHOMO_vs_DeltaELUMO_v2.mat';
    save(filename, 'el_results', 'ion_results', 'solCV_el', 'solCV_ion')
end
