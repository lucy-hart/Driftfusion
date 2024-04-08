%Use this file to symmetrically sweep HOMO/LUMO offsets vs ion
%concentration for a device with doped or undoped interlayers

%There have been some choices made here
%Doping is fixed s.t. Fermi Level offset is alway 0.1 eV from the relevant
%band edge
%Work function of the electrodes is handled as follows
%Either the Vbi is 1.1 eV, with an equal offset from the ETL and HTL CB/VB,
%respectively
%Or, if the offset is so large that a Vbi of 1.1 eV would mean the electode
%work function lying within 0.1 eV of the CB/VB edge, the work functions
%are redefined so that they lie 0.1 eV from the CB/VB edge 

%TURN SAVE OFF TO START OFF WITH (final cell)
%% Define parameter space
%Choose to use doped or undoped TLs
doped = 0;
n_values = 11;
Delta_TL = linspace(0, 0.5, n_values);
Delta_TL(1) = 1e-3;
Ion_Conc = [5e15 1e16 5e16 1e17 5e17 1e18 0];
n_ion_concs = length(Ion_Conc);

include_ions = 0;

%Rows are the values of v_sr
%Columns are the TL Energetic Offsets
params = cell(n_ion_concs, n_values);

for i=1:n_ion_concs
    for j=1:n_values
            params{i,j} = [Ion_Conc(i), Delta_TL(j)];
    end
end

%% Set up structures for storing the results
error_log = zeros(n_ion_concs, n_values);
soleq = cell(n_ion_concs, n_values);

solCV= cell(n_ion_concs, n_values);
results = cell(n_ion_concs, n_values);

%% Do (many) JV sweeps
%Select the correct input file for doped or undoped cases 
par=pc('Input_files/EnergyOffsetSweepParameters_v5_undoped_SAM.csv');


%Set the illumination for the JV sweeps 
illumination = 1;

%Reset the electrode work functions in each loop to be safe as they are
%changed for the cases where E_LUMO (E_HOMO) is far below (above) the CB
%(VB)
for i = 1:n_ion_concs
    for j = 1:n_values
        disp(["DeltaE_HOMO = ", num2str(Delta_TL(j)), " eV"])
        disp(["N_{ion} = ", num2str(Ion_Conc(i)), " cm^{-3}"])

        %HTL Energetics
        par.Phi_left = par.Phi_IP(1) + params{i,j}(2);

        %ion conc
        if i ~= n_ion_concs
            par.Ncat(:) = params{i,j}(1);
            par.Nani(:) = params{i,j}(1);
        end
        
        par = refresh_device(par);

        soleq{i,j} = equilibrate(par);
        
        %electron only scan
        if i == n_ion_concs
            Fermi_offset = par.EF0(3) - par.Phi_left;
            if Fermi_offset > 1.2
                Voc_max = Fermi_offset + 0.05;
            else
                Voc_max = 1.2;
            end
            num_points = (2*100*(Voc_max+0.2))+1;
            while Voc_max >= 1.05
                try            
                    solCV{i,j} = doCV(soleq{i, j}.el, illumination, -0.2, Voc_max, -0.2, 1e-4, 1, num_points);           
                    error_log(i,j) = 0;
                    results{i,j} = CVstats(solCV{i,j});
                    Voc_max = 0;                
                catch
                    if Voc_max > 1.05
                        warning("Electronic-only JV solution failed, reducing Vmax by 0.03 V")
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
            Fermi_offset = par.EF0(3) - par.Phi_left;
            if Fermi_offset > 1.2
                Voc_max = Fermi_offset + 0.05;
            else
                Voc_max = 1.2;
            end
            num_points = (2*100*(Voc_max+0.2))+1; 
            while Voc_max >= 1.05
                try
                    solCV{i,j} = doCV(soleq{i, j}.ion, illumination, -0.2, Voc_max, -0.2, 1e-4, 1, num_points);
                    error_log(i,j) = 0;
                    results{i,j} = CVstats(solCV{i, j});
                    Voc_max = 0;
                catch
                    if Voc_max > 1.05
                        warning("Ionic JV solution failed, reducing Vmax by 0.03 V")
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

%% Plot results 
Stats_array = zeros(n_ion_concs, n_values, 4);
e = solCV{1,1}.par.e;

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
figure('Name', 'JV Parameter vs Energy Offsets vs Ion Conc', 'Position', [50 50 800 700])
num = 2;
labels = ["J_{SC} (mA cm^{-2})", "V_{OC} (V)", "FF", "PCE (%)"];

box on 

%Coutourf has column and row inidices as x and y respectively
contourf(Delta_TL, log10(Ion_Conc(1:end-1)), Stats_array(1:end-1,:,num), 'LineWidth', 0.1)
ylabel('log_{10}(Ion Density/ cm^{-3})')
set(gca, 'YScale', 'log')
xlabel('\DeltaE_{SAM} (eV)')
xlim([0, 0.5])
%ylim([15, 18])
c = colorbar;
c.Label.String = labels(num);
%clim(lims(num,:))

%%
figure('Name', 'JV Parameter vs Energy Offsets vs Ion Conc', 'Position', [50 50 800 800])
Colours = parula(n_ion_concs-1);
num = 2;
labels = ["J_{SC} (mA cm^{-2})", "V_{OC} (V)", "FF", "PCE (%)"];
LegendLoc = ["northeast", "southwest", "southeast", "northeast"];
lims = [[-22 -10]; [0.9 1.195]; [0.15, 0.85]; [12 19]];
box on 
for i = 1:n_ion_concs
    hold on
    if i == n_ion_concs
        plot(Delta_TL, Stats_array(n_ion_concs,:,num), 'marker', 'none', 'Color', 'black')
    else
        plot(Delta_TL, Stats_array(i,:,num), 'marker', 'none', 'Color', Colours(i,:))
    end
end
xline(0.2, 'black', 'LineWidth', 2, 'LineStyle', '--')
set(gca, 'Fontsize', 25)
xlabel('Transport Layer Energetic Offset (eV)', 'FontSize', 30)
ylabel(labels(num), 'FontSize', 30)
xlim([0, 0.5])
xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
xticklabels({'0.00', '0.10', '0.20', '0.30', '0.40', '0.50'})
ylim(lims(num,:))
legend({'5e15', '1e16', '5e16', '1e17', '5e17', '1e18', 'No Ions', ''}, 'Location', LegendLoc(num), 'FontSize', 25, 'NumColumns', 2)
%title(legend, 'Ion Concentration (cm^{-3})', 'FontSize', 25)

%% Save results and solutions
save_file = 0;
if save_file == 1
    filename = 'DeltaE_v5_undoped_SAM_remix.mat'; 
    save(filename, 'results', 'solCV')
end
