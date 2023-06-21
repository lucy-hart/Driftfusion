%Use this file to symmetrically sweep HOMO/LUMO offsets vs ion
%concentration

%TURN SAVE OFF TO START OFF WITH (final cell)

tic
%% Define parameter space
%Rows are the Ion Concentrations
%Columns are the TL Energetic Offsets

n_values = 10;
Delta_TL = linspace(-0.15, 0.3, n_values);
Ion_Conc = [1e15 5e15 1e16 5e16 1e17 5e17 1e18 0];
n_ion_concs = length(Ion_Conc);

params = cell(n_ion_concs, n_values);

for i=1:n_ion_concs
    for j=1:n_values
            params{i,j} = [Ion_Conc(i), Delta_TL(j)];
    end
end

%%
error_log = zeros(n_ion_concs, n_values);
soleq = cell(n_ion_concs, n_values);
solCV = cell(n_ion_concs, n_values);
results = cell(n_ion_concs, n_values);

%% Do (many) JV sweeps
%Remeber to update the work function values if you change these parameters
%between files 
par=pc('Input_files/EnergyOffsetSweepParameters_v3.csv');
illumination = 1.1;

%Reset the electrode work functions in each loop to be safe as they are
%changed for the cases where E_LUMO (E_HOMO) is far below (above) the CB
%(VB)
for i = 1:n_ion_concs
    for j = 1:n_values
        disp(["DeltaE_HOMO = ", num2str(Delta_TL(j)), " eV"])
        if i < n_ion_concs 
            disp(["Ion_Conc = ", num2str(Ion_Conc(i)), " cm{^-3}"])
        elseif i == n_ion_concs 
            disp("No Mobile Ions")
        end
        %HTL Energetics
        par.Phi_left = -5.15;
        par.Phi_IP(1) = par.Phi_IP(3) + params{i,j}(2);
        par.Phi_EA(1) = par.Phi_IP(1) + 2.5;
        par.EF0(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
        par.Et(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
        if par.Phi_left < par.Phi_IP(1)
            par.Phi_left = par.Phi_IP(1);
        end
        %ETL Energetics
        %Need to use opposite sign at ETL to keep energy offsets symmetric
        par.Phi_right = -4.05;
        par.Phi_EA(5) = par.Phi_EA(3) - params{i,j}(2);
        par.Phi_IP(5) = par.Phi_EA(5) - 2.5;
        par.EF0(5) = (par.Phi_IP(5)+par.Phi_EA(5))/2;
        par.Et(5) = (par.Phi_IP(5)+par.Phi_EA(5))/2;
        if par.Phi_right > par.Phi_EA(5)
            par.Phi_right = par.Phi_EA(5);
        end
        %ion conc
        if i ~= n_ion_concs
            par.Ncat(:) = params{i,j}(1);
            par.Nani(:) = params{i,j}(1);
        end 

        par = refresh_device(par);

        soleq{i,j} = equilibrate(par);
        
        %electron only scan
        if i == n_ion_concs 
            Voc_max = 1.2;
            num_points = 281;
            while Voc_max >= 1.05
                try            
                    solCV{i, j} = doCV(soleq{i, j}.el, illumination, -0.2, Voc_max, -0.2, 1e-4, 1, num_points);           
                    error_log(i,j) = 0;
                    results{i,j} = CVstats(solCV{i, j});
                    Voc_max = 0;                
                catch
                    if Voc_max > 1.05
                        warning("Electronic-only JV solution failed, reducing Vmax by 0.05 V")
                        Voc_max = Voc_max - 0.03;
                        num_points = num_points - 6;
                    elseif Voc_max == 1.05
                        warning("Electronic-only JV solution failed.")
                        error_log(i,j) = 1;
                        results{i,j} = 0;
                    end
                end
            end
        
        else
            Voc_max = 1.2;
            num_points = 281; 
            while Voc_max >= 1.05
                try
                    solCV{i, j} = doCV(soleq{i, j}.ion, illumination, -0.2, Voc_max, -0.2, 1e-4, 1, num_points);
                    error_log(i,j) = 0;
                    results{i,j} = CVstats(solCV{i, j});
                    Voc_max = 0;
                catch
                    if Voc_max > 1.05
                        warning("Ionic JV solution failed, reducing Vmax by 0.05 V")
                        Voc_max = Voc_max - 0.03;
                        num_points = num_points - 6;
                    elseif Voc_max == 1.05
                        warning("Ionic JV solution failed.")
                        error_log(i,j) = 1;
                        results{i,j} = 0;
                    end
                end
            end
        end
    end
end

toc

%% Plot results 
Stats_array = zeros(n_ion_concs, n_values, 4);
for i = 1:n_ion_concs
    for j = 1:n_values
        try
            Stats_array(i,j,1) = 1e3*results{i,j}.Jsc_f;
            Stats_array(i,j,2) = results{i,j}.Voc_f;
            Stats_array(i,j,3) = results{i,j}.FF_f;
            Stats_array(i,j,4) = results{i,j}.efficiency_f;
        catch
            warning('No Stats')
            Stats_array(i,j,:) = 0;
        end
    end
end 

%%
figure('Name', 'PCE vs Energy Offsets vs Ion Conc', 'Position', [50 50 1000 1200])
Colours = parula(n_ion_concs-1);
num = 1;
labels = ["J_{SC} (mA cm^{-2})", "V_{OC} (V)", "FF", "PCE (%)"];
LegendLoc = ["northeast", "southwest", "southeast", "southeast"];
lims = [[-24 -15]; [0.77 1.24]; [0.1, 0.85]; [1 21]];
box on 
for i = 1:n_ion_concs
    hold on
    if i == 1
        plot(Delta_TL, Stats_array(n_ion_concs,:,num), 'marker', 'x', 'Color', 'black')
    else
        plot(Delta_TL, Stats_array(i-1,:,num), 'marker', 'x', 'Color', Colours(i-1,:))
    end
end
set(gca, 'Fontsize', 25)
xlabel('Transport Layer Energetic Offset (eV)', 'FontSize', 30)
ylabel(labels(num), 'FontSize', 30)
xlim([-0.15, 0.3])
ylim(lims(num,:))
legend({'No Ions', '1e15', '5e15', '1e16', '5e16', '1e17', '5e17', '1e18'}, 'Location', LegendLoc(num), 'FontSize', 25, 'NumColumns', 2)
title(legend, 'Ion Concentration (cm^{-3})', 'FontSize', 25)

%% Save results and solutions
save_file = 0;
if save_file == 1
    filename = 'tld_symmetric_DeltaEHOMO_vs_DeltaELUMO_v2.mat';
    save(filename, 'el_results', 'results', 'solCV_el', 'solCV')
end
