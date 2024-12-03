%%Load the device
read_in_par = 1;
doped = 0;

if read_in_par == 1
    if doped == 1
        par=pc('Input_files/EnergyOffsetSweepParameters_v5_doped.csv');
    elseif doped == 0
        par=pc('Input_files/EnergyOffsetSweepParameters_v5_undoped.csv');
    end
    
    Fiddle_with_Energetics = 1;
    Fiddle_with_IonConc = 1;
    IonConc = 1e18;
    %%
    if Fiddle_with_Energetics == 1
    
        DHOMO = 0.25;
        DLUMO = -0.25;
            if doped == 0
                %HTL Energetics
                par.Phi_left = -5.15;
                par.Phi_IP(1) = par.Phi_IP(3) + DHOMO;
                par.Phi_EA(1) = par.Phi_IP(1) + 2.5;
                par.EF0(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
                par.Et(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
                if par.Phi_left < par.Phi_IP(1) + 0.1
                    par.Phi_left = par.Phi_IP(1) + 0.1;
                end
    
                %ETL Energetics
                par.Phi_right = -4.05;
                par.Phi_EA(5) = par.Phi_EA(3) + DLUMO;
                par.Phi_IP(5) = par.Phi_EA(5) - 2.5;
                par.EF0(5) = (par.Phi_IP(5)+par.Phi_EA(5))/2;
                par.Et(5) = (par.Phi_IP(5)+par.Phi_EA(5))/2;
                if par.Phi_right > par.Phi_EA(5) - 0.1
                    par.Phi_right = par.Phi_EA(5) - 0.1;
                end
            
            elseif doped == 1
                %HTL Energetics
                par.Phi_left = -5.15;
                par.Phi_IP(1) = par.Phi_IP(3) + DHOMO;
                par.Phi_EA(1) = par.Phi_IP(1) + 2.5;
                par.EF0(1) = par.Phi_IP(1) + 0.1;
                par.Et(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
                if par.Phi_left < par.Phi_IP(1) + 0.1
                    par.Phi_left = par.Phi_IP(1) + 0.1;
                end
    
                %ETL Energetics
                %Need to use opposite sign at ETL to keep energy offsets symmetric
                par.Phi_right = -4.05;
                par.Phi_EA(5) = par.Phi_EA(3) + DLUMO;
                par.Phi_IP(5) = par.Phi_EA(5) - 2.5;
                par.EF0(5) = par.Phi_EA(5) - 0.1;
                par.Et(5) = (par.Phi_IP(5) + par.Phi_EA(5))/2;
                if par.Phi_right > par.Phi_EA(5) - 0.1
                    par.Phi_right = par.Phi_EA(5) - 0.1;
                end
    
            end
        
        par = refresh_device(par);
    
    end
    
    if Fiddle_with_IonConc == 1
    
       par.Ncat(:) = IonConc;
       par.Nani(:) = IonConc;
    
       par = refresh_device(par);
    
    end
    
    par.vsr_mode = 1;
    par.frac_vsr_zone = 0.05;
    par = refresh_device(par);
    eqm = equilibrate(par);
end
%% Get capacitance 
light_intensity = 1;
%NB: haven't use exactly V_flat as devices don't stabilise here for some reason.
if doped == 1
    %Vbias_ar = [0 0.9001 1.12];
    %Vbias_ar = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.901 1.0 1.1];
elseif doped == 0
    %Vbias_ar = [0];
    Vbias_ar = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.901 1.0 1.1];
end
ions_ar = [0 1];
tstab = 120;
DeltaV = 0.01;
q = 1.6e-19;

results = zeros(numel(Vbias_ar), 7);

%% Stabilise device at initial voltage (V)
for i = 1:numel(ions_ar)
    ions = ions_ar(i);
    for j = 1:numel(Vbias_ar)
        Vbias = Vbias_ar(j);
        if light_intensity ~= 0
            if ions == 1
                sol_ill = changeLight(eqm.ion, light_intensity, 0, 1);
            elseif ions == 0
                sol_ill = changeLight(eqm.el, light_intensity, 0, 1);
            end
        else
            if ions == 1
                sol_ill = eqm.ion;
            elseif ions == 0
                sol_ill = eqm.el;
            end
        end
        
        par = sol_ill.par;
        
        %ramp voltage up to the applied voltage 
        par.tmesh_type = 1;
        par.t0 = 0;
        par.tmax = 1e-2;
        par.tpoints = 100;
        
        par.V_fun_type = 'sweep';
        par.V_fun_arg(1) = 0;
        par.V_fun_arg(2) = Vbias;
        par.V_fun_arg(3) = 1e-2;
        
        sol = df(sol_ill, par);
        
        par = sol.par;
        
        %Hold device at Vbias for time tstab
        par.tmesh_type = 1;
        par.t0 = 0;
        par.tmax = tstab;
        par.tpoints = 100;
        
        par.V_fun_type = 'constant';
        par.V_fun_arg(1) = Vbias;
            
        sol = df(sol, par);
        
        %% Calculate the charge per unit area (Q m-2)
        x = sol.x*1e-2;
        
        n_el = 1e6*(sol.u(end,:,3) - sol.u(end,:,2));
        n_bar = 1e6*(sol.u(end,:,3) + sol.u(end,:,2)).*0.5;
        
        if ions == 1
            Q_el_electrode1 = q*abs(trapz(x, n_el));
            Q_el1 = q*trapz(x, abs(n_el));
            Q_chem1 = q*trapz(x, n_bar);
            n_ion = 1e6*(sol.u(end,:,4) - sol.u(end,1,4));
            Q_ion1 = q*trapz(x, abs(n_ion));
            Q_net1 = q*trapz(x, abs(n_ion+n_el));
        elseif ions == 0
            Q_el_electrode_noions1 = q*abs(trapz(x, n_el));
            Q_el_noions1 = q*trapz(x, abs(n_el));
            Q_chem_noions1 = q*trapz(x, n_bar);
        end
        
        %% Stabilise device at new voltage (V+DeltaV)
        par = sol_ill.par;
        
        %ramp voltage up to the applied voltage 
        par.tmesh_type = 1;
        par.t0 = 0;
        par.tmax = 1e-2;
        par.tpoints = 100;
        
        par.V_fun_type = 'sweep';
        par.V_fun_arg(1) = 0;
        par.V_fun_arg(2) = Vbias+DeltaV;
        par.V_fun_arg(3) = 1e-2;
        
        sol = df(sol_ill, par);
        
        par = sol.par;
        
        %Hold device at new Vbias for time tstab
        par.tmesh_type = 1;
        par.t0 = 0;
        par.tmax = tstab;
        par.tpoints = 100;
        
        par.V_fun_type = 'constant';
        par.V_fun_arg(1) = Vbias+DeltaV;
            
        sol = df(sol, par);
        
        %% Calculate the new charge per unit area (Q m-2)
        x = sol.x*1e-2;
        
        n_el = 1e6*(sol.u(end,:,3) - sol.u(end,:,2));
        n_bar = 1e6*(sol.u(end,:,3)+sol.u(end,:,2)).*0.5;
        
        if ions == 1
            Q_el_electrode2 = q*abs(trapz(x, n_el));
            Q_el2 = q*trapz(x, abs(n_el));
            Q_chem2 = q*trapz(x, n_bar);
            n_ion = 1e6*(sol.u(end,:,4) - sol.u(end,1,4));
            Q_ion2 = q*trapz(x, abs(n_ion));
            Q_net2 = q*trapz(x, abs(n_ion+n_el));
        elseif ions == 0
            Q_el_electrode_noions2 = q*abs(trapz(x, n_el));
            Q_el_noions2 = q*trapz(x, abs(n_el));
            Q_chem_noions2 = q*trapz(x, n_bar);
        end
        
        %% Get the capacitance 
        %Have converted distances so units are F m-2
        if ions == 1
            results(j,1) = abs(Q_el_electrode2 - Q_el_electrode1)/DeltaV;
            results(j,2) = abs(Q_el2 - Q_el1)/DeltaV;
            results(j,3) = abs(Q_net2 - Q_net1)/DeltaV;
            results(j,4) = abs(Q_chem2 - Q_chem1)/DeltaV;
        elseif ions == 0          
            results(j,5) = abs(Q_el_electrode_noions2 - Q_el_electrode_noions1)/DeltaV;
            results(j,6) = abs(Q_el_noions2 - Q_el_noions1)/DeltaV;
            results(j,7) = abs(Q_chem_noions2 - Q_chem_noions1)/DeltaV;
        end

    end
end

