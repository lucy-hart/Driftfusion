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

%Used this file to produce v4 of the simulation results 

%TURN SAVE OFF TO START OFF WITH (final cell)

%% Define parameter space
%Choose to use doped or undoped TLs
doped = 0;
n_values = 7;
Delta_TL = linspace(0, 0.3, n_values);
Symmetric_offset = 0;
%Fix the offset for the ETL or HTL
Fix_ETL = 1;
%Energetic offset between the perovskite and TL for the TL with fixed
%energetics 
Fixed_offset = 1e-5;
%This is a bit of a hack, but if the offfset is exactly 0, the surface
%recombination error becomes huge for reasons I do not fully understand...
%The saga continues - only seems to matter for the HTL, not the ETL...
Delta_TL(1) = 1e-5;
Ion_Conc = 1e17;

%% Set up structures for storing the results
error_log = zeros(2,n_values);
soleq = cell(2,n_values);
solCV = cell(2,n_values);
results = cell(2,n_values);

%% Do (many) JV sweeps
%Select the correct input file for doped or undoped cases 
if doped == 1
    par=pc('Input_files/EnergyOffsetSweepParameters_v5_doped.csv');
elseif doped == 0
    par=pc('Input_files/EnergyOffsetSweepParameters_v5_undoped.csv');
end 

%Set the illumination for the JV sweeps 
illumination = 1;

%Set Ion Concentration
par.Ncat(:) = Ion_Conc;
par.Nani(:) = Ion_Conc;

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
for i = 1:2
    for j = 1:n_values
        disp(["DeltaE_HOMO = ", num2str(Delta_TL(j)), " eV"])
        if i == 1 
            disp(["Ion_Conc = ", num2str(Ion_Conc), " cm{^-3}"])
        elseif i == 2 
            disp("No Mobile Ions")
        end

        %HTL Energetics
        if Symmetric_offset == 1 || (Symmetric_offset == 0 && Fix_ETL == 1)
            par.Phi_left = -5.15;
            par.Phi_IP(1) = par.Phi_IP(3) + Delta_TL(j);
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
        elseif Symmetric_offset == 0 && Fix_ETL == 0
            par.Phi_left = -5.15;
            par.Phi_IP(1) = par.Phi_IP(3) + Fixed_offset;
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
        end

        %ETL Energetics
        %Need to use opposite sign at ETL to keep energy offsets symmetric
        if Symmetric_offset == 1 || (Symmetric_offset == 0 && Fix_ETL == 0)
            par.Phi_right = -4.05;
            par.Phi_EA(5) = par.Phi_EA(3) - Delta_TL(j);
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
        elseif Symmetric_offset == 0 && Fix_ETL == 1
            par.Phi_right = -4.05;
            par.Phi_EA(5) = par.Phi_EA(3) - Fixed_offset;
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
        end
        
        par = refresh_device(par);

        soleq{i,j} = equilibrate(par);
        
        %electron only scan
        if i == 2
            Fermi_offset = par.EF0(5) - par.EF0(1);
            if Fermi_offset > 1.2
                Voc_max = Fermi_offset + 0.05;
            else
                Voc_max = 1.2;
            end
            num_points = (2*100*(Voc_max+0.2))+1;
            while Voc_max >= 1.05
                try            
                    solCV{i, j} = doCV(soleq{i, j}.el, illumination, -0.2, Voc_max, -0.2, 1e-4, 1, num_points);           
                    error_log(i,j) = 0;
                    results{i,j} = CVstats(solCV{i, j});
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
            Fermi_offset = par.EF0(5) - par.EF0(1);
            if Fermi_offset > 1.2
                Voc_max = Fermi_offset + 0.05;
            else
                Voc_max = 1.2;
            end
            num_points = (2*100*(Voc_max+0.2))+1; 
            while Voc_max >= 1.05
                try
                    solCV{i, j} = doCV(soleq{i, j}.ion, illumination, -0.2, Voc_max, -0.2, 1e-4, 1, num_points);
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

%%
%Find V_invert for the ionic solutions
V_invert = zeros(1,n_values);

for i = 1:n_values
    temp_value = findVinvert(solCV{1,i});
    temp_val_HTL = temp_value.HTL;
    temp_val_ETL = temp_value.ETL;
    if temp_val_ETL ~= 0 && temp_val_HTL ~=0
        V_invert(1,i) = min(temp_val_HTL, temp_val_ETL);
    elseif temp_val_ETL == 0
        V_invert(1,i) = temp_val_HTL;
    elseif temp_val_HTL == 0
        V_invert(1,i) = temp_val_ETL;
    end 
end

%%
%Do fixed ion JV at V_pre = V_invert_min
fixed_ion_JVs = cell(1,n_values);
J_fixed_ion = cell(1,n_values);
J_el = cell(1,n_values);

tstab = 120;

for i = 1:n_values
    sol_ill = changeLight(soleq{1,i}.ion, illumination, 0, 1);
    par = sol_ill.par;
    par.tmesh_type = 1;
    par.t0 = 0;
    par.tmax = 1e-2;
    par.tpoints = 100;

    par.V_fun_type = 'sweep';
    par.V_fun_arg(1) = 0;
    par.V_fun_arg(2) = V_invert(i);
    par.V_fun_arg(3) = 1e-2;

    sol = df(sol_ill, par);

    par = sol.par;

    %Hold device at Vbias for time tstab
    par.tmesh_type = 1;
    par.t0 = 0;
    par.tmax = tstab;
    par.tpoints = 100;

    par.V_fun_type = 'constant';
    par.V_fun_arg(1) = V_invert(i);
    
    disp(['Stabilising solution at ' num2str(V_invert(i)) ' V'])
    sol = df(sol, par);

    disp(['Doing JV for Vstab = ' num2str(V_invert(i)) ' V'])
    sol.par.mobseti = 0;
    fixed_ion_JVs{i} = doCV(sol, illumination, -0.2, 1.2, -0.2, 1e-3, 1, 281);
    J_fixed_ion{i} = dfana.calcJ(fixed_ion_JVs{i}).tot(:,1);
    J_el{i} = dfana.calcJ(solCV{2,i}).tot(:,1);
    if i == 1
            V_fixed_ion = dfana.calcVapp(fixed_ion_JVs{i});
    end
end

%% Plot electron only JVs and fixed ion JVs at V_bias = min(V_invert)
figure('Name', 'Asymmetric V_invert Check')
cmap = colormap(parula(length(V_invert)));
labels = {'0.00', '0.05', '0.10', '0.15', '0.20', '0.25', '0.30'};
hold on
box on
xline(0, 'black', 'HandleVisibility', 'off')
yline(0, 'black', 'HandleVisibility', 'off')
for i = 1:length(V_invert)
    plot(V_fixed_ion, 1e3*J_fixed_ion{i}, 'DisplayName', labels{i}, 'color', cmap(i,:), 'LineStyle', '-')
    V_el = dfana.calcVapp(solCV{2,i});
    num_points = length(V_el);
    stop = int32(ceil(num_points/2));
    plot(V_el(1:stop), 1e3*J_el{i}(1:stop), 'color', 'black', 'LineStyle', '--', 'HandleVisibility', 'Off')
end
hold off
ylabel('Current Density (mA cm^{-2})')
ylim([-25, 10])
xlabel('Voltage (V)')
xlim([0, 1.2])
legend()
title(legend, '\Delta E_{TL} (eV)')

%%
V_bi = zeros(1, n_values);
DeltaEF_TL = zeros(1, n_values);

for i = 1:n_values
    V_bi(i) = soleq{1,i}.ion.par.Phi_right - soleq{1,i}.ion.par.Phi_left;
    DeltaEF_TL(i) = soleq{1,i}.ion.par.Phi_EA(5) - soleq{1,i}.ion.par.Phi_IP(1);
end