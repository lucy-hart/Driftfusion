function J_return = doSaP_multi_v2(sol_ini, tstab, Vbias, Vpulse, light_intensity, verbose)
% Performs a simulation of a stabilise and pulse (SaP) measurement
% Input arguments:
% SOL_INI = solution containing intitial conditions (dark eqm device)
% VBIAS = array of stabilisation biases to sample at
% VPULSE = array of voltages to sample at in the pulsed JV
% TRAMP = the rise time of the voltage pulse (8e-4 s from paper)
% TSAMPLE = time after voltage pulse when current is measured
% TSTAB = time for which the device is stabilised at Vbias
% LIGHT_INTENSITY = light intensity at which the measurement is done

%% Start Code
disp('Starting SaP')

num_bias = length(Vbias);
num_pulses = length(Vpulse);

%Each row in the solution structure is a different prebais
%The first column is the solution after being held at the relevant Vbias
%for time tstab
%The other columns are the resuls of the pulsed JVs
J_return = zeros(num_bias, num_pulses);

%% Go to correct light intensity
%Funtion assumes you are using the first light source currently 
sol_ill = changeLight(sol_ini, light_intensity, 0, 1);

%% 
for i = 1:num_bias    
    %Generate solutions after voltage stabilisation period
    par = sol_ill.par;
    par.vsr_check = 0;
    
    try
        %ramp voltage up to the applied voltage 
        par.tmesh_type = 1;
        par.t0 = 0;
        par.tmax = 1e-2;
        par.tpoints = 100;
    
        par.V_fun_type = 'sweep';
        par.V_fun_arg(1) = 0;
        par.V_fun_arg(2) = Vbias(i);
        par.V_fun_arg(3) = 1e-2;
    
        sol = df(sol_ill, par);
    
        par = sol.par;
        par.vsr_check = 0;
    
        %Hold device at Vbias for time tstab
        par.tmesh_type = 1;
        par.t0 = 0;
        par.tmax = tstab;
        par.tpoints = 100;
    
        par.V_fun_type = 'constant';
        par.V_fun_arg(1) = Vbias(i);
        
        if verbose == 1
            disp(['Stabilising solution at ' num2str(Vbias(i)) ' V'])
        end
        SaPsol = df(sol, par);
    catch 
        try
            %ramp voltage up to the applied voltage 
            par.tmesh_type = 1;
            par.t0 = 0;
            par.tmax = 1e-2;
            par.tpoints = 100;
        
            par.V_fun_type = 'sweep';
            par.V_fun_arg(1) = 0;
            par.V_fun_arg(2) = Vbias(i) + 0.01;
            par.V_fun_arg(3) = 1e-2;
        
            sol = df(sol_ill, par);
        
            par = sol.par;
            par.vsr_check = 0;
        
            %Hold device at Vbias for time tstab
            par.tmesh_type = 1;
            par.t0 = 0;
            par.tmax = tstab;
            par.tpoints = 100;
        
            par.V_fun_type = 'constant';
            par.V_fun_arg(1) = Vbias(i) + 0.01;
            
            if verbose == 1
                disp(['Stabilising solution at ' num2str(Vbias(i) + 0.01) ' V'])
            end
            SaPsol = df(sol, par);
    
            Vbias(i) = Vbias(i) + 0.01;
        catch
            try
                %ramp voltage up to the applied voltage 
                par.tmesh_type = 1;
                par.t0 = 0;
                par.tmax = 1e-2;
                par.tpoints = 100;
            
                par.V_fun_type = 'sweep';
                par.V_fun_arg(1) = 0;
                par.V_fun_arg(2) = Vbias(i) - 0.01;
                par.V_fun_arg(3) = 1e-2;
            
                sol = df(sol_ill, par);
            
                par = sol.par;
                par.vsr_check = 0;
            
                %Hold device at Vbias for time tstab
                par.tmesh_type = 1;
                par.t0 = 0;
                par.tmax = tstab;
                par.tpoints = 100;
            
                par.V_fun_type = 'constant';
                par.V_fun_arg(1) = Vbias(i) - 0.01;
                
                if verbose == 1
                    disp(['Stabilising solution at ' num2str(Vbias(i) - 0.01) ' V'])
                end
                SaPsol = df(sol, par);
        
                Vbias(i) = Vbias(i) - 0.01;
            catch
            warning(['Could not stabilise device at ' num2str(Vbias(i)) ' V'])
            Vbias(i) = 100;
            end
        end 
    end

    if Vbias(i) ~= 100
        SaPsol.par.vsr_check = 0;
        %turn off ion motion for the duration of the JV
        %assuming that this is a valid assumption
        SaPsol.par.mobseti = 0;
        try 
            fixed_ion_JV = doCV(SaPsol, light_intensity, -0.2, 1.2, -0.2, 1e-2, 0.5, 146);
            V_fixed_ion = dfana.calcV(fixed_ion_JV);
            J_fixed_ion = dfana.calcJ(fixed_ion_JV).tot(:,1);
            J_return(i,:) = interp1(V_fixed_ion(1:end), J_fixed_ion(1:end), Vpulse);
        catch 
            warning(['Could not do JV for Vbias = ' num2str(Vbias(i)) ' V'])
            J_return(i,:) = 0;
        end 
    else
        J_return(i,:) = 0;
    end

end    

