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
doped = 0;
high_performance = 0;
tau_n = [1e-6 1e-7 1e-8 1e-9];
E_trap = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7];
n_trap  = length(E_trap);
n_values = length(tau_n);
%Rows are the values of E_trap
%Columns are the TL Energetic Offsets
params = cell(n_trap, n_values);

for i=1:n_trap
    for j=1:n_values
            params{i,j} = [E_trap(i), tau_n(j)];
    end
end

%% Set up structures for storing the results
error_log = zeros(n_trap, n_values);
soleq = cell(n_trap, n_values);

solCV_ion = cell(n_values, n_values);
results_ion = cell(n_values, n_values);

solCV_el = cell(n_values, n_values);
results_el = cell(n_values, n_values);

%% Do (many) JV sweeps
%Select the correct input file for doped or undoped cases 
if doped == 1
    if high_performance == 0
        par=pc('Input_files/EnergyOffsetSweepParameters_v5_doped.csv');
    elseif high_performance == 1
        par=pc('Input_files/EnergyOffsetSweepParameters_v5_doped_higherPCE.csv');
    end
elseif doped == 0
    if high_performance == 0
        par=pc('Input_files/EnergyOffsetSweepParameters_v5_undoped_Weidong.csv');
    elseif high_performance == 1
        par=pc('Input_files/EnergyOffsetSweepParameters_v5_undoped_higherPCE.csv');
    end
end

%Set the illumination for the JV sweeps and the max voltage to go up to in
%the JV sweeps
if high_performance == 0
    illumination = 1;
    max_val = 1.2;
    phi_L = -5.15;
    phi_R = -4.05;
elseif high_performance == 1
    illumination = 1.2;
    max_val = 1.3;
    phi_L = -5.2;
    phi_R = -4.1;
end

%Reset the electrode work functions in each loop to be safe as they are
%changed for the cases where E_LUMO (E_HOMO) is far below (above) the CB
%(VB)
for i = 1:n_trap
    for j = 1:n_values
        disp(["tau_n = ", num2str(tau_n(j)), " s"])
        disp(["Trap depth = ", num2str(E_trap(i)), " eV"])

%         par.Et(2:4) = (par.Phi_IP(3) + par.Phi_EA(3))/2 + params{i,j}(1);
        par.Et(2:4) = (par.Phi_IP(3) + par.Phi_EA(3))/2 - params{i,j}(1);
        par.taun(3) = params{i,j}(2);
        par.taup(3) = 100e-9;

        %ion conc
        par.Ncat(:) = 1e18;
        par.Nani(:) = 1e18;
        
        %Do this as it seesm to reduce the discrepency between surface
        %volumetric surace recombination model and the abrupt Thinterface one
        %But also made solution less stable? Maybe better to tinker with
        %this on the one which varies surface recombination...
        %par.frac_vsr_zone = 0.05;
        par = refresh_device(par);

        soleq{i,j} = equilibrate(par);
        
        %electron only scan
        Fermi_offset = par.EF0(5) - par.EF0(1);
        if Fermi_offset > 1.2
            Voc_max = Fermi_offset + 0.05;
        else
            Voc_max = 1.2;
        end
        num_points = (2*100*(Voc_max+0.2))+1;
        while Voc_max >= 1.05
            try            
                solCV_el{i,j} = doCV(soleq{i, j}.el, illumination, -0.2, Voc_max, -0.2, 1e-4, 1, num_points);           
                error_log(i,j) = 0;
                results_el{i,j} = CVstats(solCV_el{i,j});
                Voc_max = 0;                
            catch
                if Voc_max > 1.05
                    warning("Electronic-only JV solution failed, reducing Vmax by 0.03 V")
                    Voc_max = Voc_max - 0.03;
                    num_points = num_points - 6;
                elseif Voc_max == 1.05
                    warning("Electronic-only JV solution failed.")
                    error_log(i,j) = 1;
                    results_el{i,j} = 0;
                end
            end
        end
    
        Fermi_offset = par.EF0(5) - par.EF0(1);
        if Fermi_offset > 1.2
            Voc_max = Fermi_offset + 0.05;
        else
            Voc_max = 1.2;
        end
        num_points = (2*100*(Voc_max+0.2))+1; 
        while Voc_max >= 1.05
            try
                solCV_ion{i,j} = doCV(soleq{i, j}.ion, illumination, -0.2, Voc_max, -0.2, 1e-4, 1, num_points);
                error_log(i,j) = 0;
                results_ion{i,j} = CVstats(solCV_ion{i, j});
                Voc_max = 0;
            catch
                if Voc_max > 1.05
                    warning("Ionic JV solution failed, reducing Vmax by 0.03 V")
                    Voc_max = Voc_max - 0.03;
                    num_points = num_points - 6;
                elseif Voc_max == 1.05
                    warning("Ionic JV solution failed.")
                    error_log(i,j) = 1;
                    results_ion{i,j} = 0;
                end
            end
        end
    end
end

%% Plot results 
Stats_array_el = zeros(n_trap, n_values, 4);
Stats_array_ion = zeros(n_trap, n_values, 4);
e = solCV_ion{1,1}.par.e;

for i = 1:n_trap
    for j = 1:n_values
        try
            Stats_array_el(i,j,1) = 1e3*results_el{i,j}.Jsc_f;
            Stats_array_el(i,j,2) = results_el{i,j}.Voc_f;
            Stats_array_el(i,j,3) = results_el{i,j}.FF_f;
            Stats_array_el(i,j,4) = results_el{i,j}.efficiency_f;
        catch
            warning('No Stats')
            Stats_array_el(i,j,:) = 0;
        end
    end
end 

for i = 1:n_trap
    for j = 1:n_values
        try
            Stats_array_ion(i,j,1) = 1e3*results_ion{i,j}.Jsc_f;
            Stats_array_ion(i,j,2) = results_ion{i,j}.Voc_f;
            Stats_array_ion(i,j,3) = results_ion{i,j}.FF_f;
            Stats_array_ion(i,j,4) = results_ion{i,j}.efficiency_f;
        catch
            warning('No Stats')
            Stats_array_ion(i,j,:) = 0;
        end
    end
end  

%%
figure('Name', 'JV Parameter vs Energy Offsets vs Ion Conc', 'Position', [50 50 800 700])
num = 1;
labels = ["J_{SC} (mA cm^{-2})", "V_{OC} (V)", "FF", "PCE (%)"];

if doped == 1
    lims = [[-22 -18]; [0.8 1.195]; [0.15, 0.85]; [7 21]];
else
    lims = [[-22 -19]; [0.55 1.16]; [0.15, 0.85]; [5 25.5]];
end

box on 

%Coutourf has column and row inidices as x and y respectively
contourf(log10(tau_n), 0.8*ones(1,length(E_trap)) - E_trap, Stats_array_ion(:,:,num)-Stats_array_el(:,:,num), 'LineWidth', 0.1)
ylabel('Surface Recombination Velocity (cm s^{-1})')
%set(gca, 'YScale', 'log')
xlabel('\DeltaE_{TL} (eV)')
% xlim([0, 0.3])
% ylim([1, 1e3])
c = colorbar;
c.Label.String = labels(num);
%clim(lims(num,:))

%%
figure('Name', 'JV Parameter vs Energy Offsets vs Ion Conc', 'Position', [50 50 800 700])
num = 1;
labels = ["J_{SC} (mA cm^{-2})", "V_{OC} (V)", "FF", "PCE (%)"];
LegendLoc = ["northeast", "southwest", "southeast", "northeast"];

box on 

%Coutourf has column and row inidices as x and y respectively
contourf(Delta_TL, 0.8 - E_trap, Stats_array_ion(:,:,num)-Stats_array_el(:,:,num), 'LineWidth', 0.1)
ylabel('Surface Recombination Velocity (cm s^{-1})')
set(gca, 'YScale', 'log')
xlabel('\DeltaE_{TL} (eV)')
xlim([0, 0.3])
ylim([1, 1e3])
c = colorbar;
c.Label.String = labels(num);
%clim([-0.1,2])

%% Save results and solutions
save_file = 0;
if save_file == 1
    if doped == 0
        filename = 'DeltaE_v5_undoped.mat';
    elseif doped == 1
        filename = 'DeltaE_v5_doped_VaryHTL_FixedETL_0p15eV_vsr_100cms-1.mat';
    end 
    save(filename, 'results', 'solCV')
end
