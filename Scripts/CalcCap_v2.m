%%Load the device
read_in_par = 0;
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
light_intensity = 0;
Vbias = 0;
ions = 1;
tstab = 120;
DeltaV = 0.01;
q = 1.6e-19;

results = zeros(numel(Vbias_ar), 2);

%% Stabilise device at initial voltage (V)

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

%Hold device at new Vbias for time tstab
par.tmesh_type = 1;
par.t0 = 0;
par.tmax = tstab;
par.tpoints = 100;

par.V_fun_type = 'constant';
par.V_fun_arg(1) = Vbias;
    
sol = df(sol, par);

%% Apply the voltage step
par = sol.par;

%apply voltage pulse
par.tmesh_type = 1;
par.t0 = 0;
par.tmax = 100;
par.tpoints = 2000;

par.V_fun_type = 'smoothed_square';
par.V_fun_arg(1) = Vbias;
par.V_fun_arg(2) = Vbias+DeltaV;
par.V_fun_arg(3) = 150;
par.V_fun_arg(4) = 95;

sol = df(sol, par);
J = dfana.calcJ(sol).tot;
t = sol.t;


%%
J_trans = J(:,1) - J(end,1);
argstart = find(diff(J_trans>=0),1);
Q = trapz(t(argstart:end), J_trans(argstart:end)); 
Cap = 1e4*Q/DeltaV;

