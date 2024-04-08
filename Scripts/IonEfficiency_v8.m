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
n_values = 7;
Delta_TL = linspace(0, 0.3, n_values);
%This is a bit of a hack, but if the offfset is exactly 0, the surface
%recombination error becomes huge for reasons I do not fully understand...
%The saga continues - only seems to matter for the HTL, not the ETL...
Delta_TL(1) = 1e-3;

%Rows are the values of v_sr
%Columns are the TL Energetic Offsets
params = cell(n_values, n_values);
for i=1:n_values
    for j=1:n_values
            params{i,j} = [Delta_TL(i), Delta_TL(j)];
    end
end

%% Set up structures for storing the results
error_log = zeros(n_values, n_values);
soleq = cell(n_values, n_values);

solCV_ion = cell(n_values, n_values);
results_ion = cell(n_values, n_values);

solCV_el = cell(n_values, n_values);
results_el = cell(n_values, n_values);

%% Do (many) JV sweeps
%Select the correct input file for doped or undoped cases 
if doped == 1
    par=pc('Input_files/EnergyOffsetSweepParameters_v5_doped.csv');
elseif doped == 0
    par=pc('Input_files/EnergyOffsetSweepParameters_v5_undoped.csv');
end 

%Set the illumination for the JV sweeps 
illumination = 1;

%Piers version - set this to 1 to ensure that the WFs of the electrodes is
%always equal to the Fermi level of the contacts - stops there being a
%Schottky barrier at the TL/metal interface, which negatively affects
%device performance. 
%Less similar to the experimental situation maybe, as there the same
%electrode materials are used for both interlayers
Piers_version = 0;
%Only makes sense to do this for the doped case so put this is as a
%protection in case you forget to change it
if doped == 0
    Piers_version = 0;
end

%Reset the electrode work functions in each loop to be safe as they are
%changed for the cases where E_LUMO (E_HOMO) is far below (above) the CB
%(VB)
for i = 1:n_values
    for j = 1:n_values
        disp(["DeltaE_HOMO = ", num2str(Delta_TL(j)), " eV"])
        disp(["DeltaE_LUMO = ", num2str(Delta_TL(i)), " eV"])

        %HTL Energetics        
        par.Phi_left = -5.15;
        par.Phi_IP(1) = par.Phi_IP(3) + params{i,j}(2);
        par.Phi_EA(1) = par.Phi_IP(1) + 2.5;
        par.Et(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
        if doped == 0
            par.EF0(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
        elseif doped == 1
            par.EF0(1) = par.Phi_IP(1) + 0.1;
        end
        if Piers_version == 0
            if par.Phi_left < par.Phi_IP(1) + 0.1
                par.Phi_left = par.Phi_IP(1) + 0.1;
            end
        elseif Piers_version == 1
            par.Phi_left = par.Phi_IP(1) + 0.1;
        end


        %ETL Energetics
        %Need to use opposite sign at ETL to keep energy offsets symmetric
        par.Phi_right = -4.05;
        par.Phi_EA(5) = par.Phi_EA(3) - params{i,j}(1);
        par.Phi_IP(5) = par.Phi_EA(5) - 2.5;
        par.Et(5) = (par.Phi_IP(5) + par.Phi_EA(5))/2;
        if doped == 0
            par.EF0(5) = (par.Phi_IP(5) + par.Phi_EA(5))/2;
        elseif doped == 1
            par.EF0(5) = par.Phi_EA(5) - 0.1;
        end
        if Piers_version == 0
            if par.Phi_right > par.Phi_EA(5) - 0.1
                par.Phi_right = par.Phi_EA(5) - 0.1;
            end
        elseif Piers_version == 1
            par.Phi_right = par.Phi_EA(5) - 0.1;
        end


        %ion conc
        par.Ncat(:) = 1e18;
        par.Nani(:) = 1e18;
        
        %Do this as it seesm to reduce the discrepency between surface
        %volumetric surace recombination model and the abrupt interface one
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
                solCV_el{i,j} = doCV(soleq{i,j}.el, illumination, -0.2, Voc_max, -0.2, 1e-4, 1, num_points);           
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
                    error_log(i,j) = error_log(i,j) + 1;
                    results_el{i,j} = 0;
                end
            end
        end
    
        %Ionic scan
        Fermi_offset = par.EF0(5) - par.EF0(1);
        if Fermi_offset > 1.2
            Voc_max = Fermi_offset + 0.05;
        else
            Voc_max = 1.2;
        end
        num_points = (2*100*(Voc_max+0.2))+1; 
        while Voc_max >= 1.05
            try
                solCV_ion{i,j} = doCV(soleq{i,j}.ion, illumination, -0.2, Voc_max, -0.2, 1e-4, 1, num_points);
                error_log(i,j) = 0;
                results_ion{i,j} = CVstats(solCV_ion{i,j});
                Voc_max = 0;
            catch
                if Voc_max > 1.05
                    warning("Ionic JV solution failed, reducing Vmax by 0.03 V")
                    Voc_max = Voc_max - 0.03;
                    num_points = num_points - 6;
                elseif Voc_max == 1.05
                    warning("Ionic JV solution failed.")
                    error_log(i,j) = error_log(i,j) + 1;
                    results_ion{i,j} = 0;
                end
            end
        end
    end
end

%% Plot results 
Stats_array_el = zeros(n_values, n_values, 4);
Stats_array_ion = zeros(n_values, n_values, 4);
e = solCV_ion{1,1}.par.e;

for i = 1:n_values
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

for i = 1:n_values
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
num = 2;
labels = ["J_{SC} (mA cm^{-2})", "V_{OC} (V)", "FF", "PCE (%)"];

if doped == 1
    lims = [[-22 -18]; [0.98 1.165]; [0.15, 0.85]; [15 21]];
else
    lims = [[-22 -18]; [0.8 1.155]; [0.15, 0.85]; [13 19]];
end

box on 

%Coutourf has column and row inidices as x and y respectively
%Data has the form solCV{Delta_HTL, Delta_ETL}
%For doped, light comes through ETL
%For undoped "             "  HTL
contourf(Delta_TL, Delta_TL, Stats_array_el(:,:,num), 'LineWidth', 0.1)
xlabel('\DeltaE_{ETL} (eV)')
ylabel('\DeltaE_{HTL} (eV)')
xlim([0, 0.3])
ylim([0, 0.3])
c = colorbar;
c.Label.String = labels(num);
clim(lims(num,:))


%% Save results and solutions
save_file = 0;
if save_file == 1
    if doped == 0
        filename = 'DeltaE_v5_undoped.mat';
    elseif doped == 1
        filename = 'DeltaETL_BothSides_v5.mat';
    end 
    save(filename, 'results', 'solCV')
end

