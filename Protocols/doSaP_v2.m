function SaPsol = doSaP_v2(sol_ini, Vbias, Vpulse, tramp, tsample, tstab, light_intensity)
% Performs a simulation of a stabilise and pulse (SaP) measurement
% Input arguments:
% SOL_INI = solution containing intitial conditions (dark eqm device)
% VBIAS = array of stabilisation biases to sample at
% VPULSE = array of voltages to sample at in the pulsed JV
% TRAMP = the rise time of the voltage pulse (8e-4 s from paper)
% TSAMPLE = time after voltage pulse when current is measured
% TSTAB = time for which the device is stabilised at Vbias
% LIGHT_INTENSITY = light intensity at which the measurement is done
% 
%

%% Start Code
disp('Starting SaP')

num_bias = length(Vbias);
num_pulses = length(Vpulse);

%Each row in the solution structure is a different prebais
%The first column is the solution after being held at the relevant Vbias
%for time tstab
%The other columns are the resuls of the pulsed JVs
SaPsol = cell(num_bias, num_pulses+1);

%% Go to correct light intensity
%Funtion assumes you are using the first light source currently 
sol_ill = changeLight(sol_ini, light_intensity, 0, 1);

%% Generate solutions after voltage stabilisation period
for i = 1:num_bias
    
    par = sol_ill.par;

    if Vbias(i) >= 1
        par.mu_p(1) = par.mu_p(1)/1000;
    end
    
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
    
        %Hold device at Vbias for time tstab
        par.tmesh_type = 1;
        par.t0 = 0;
        par.tmax = tstab;
        par.tpoints = 100;
    
        par.V_fun_type = 'constant';
        par.V_fun_arg(1) = Vbias(i);
        
        disp(['Stabilising solution at ' num2str(Vbias(i)) ' V'])

        SaPsol{i,1} = df(sol, par);
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
        
            %Hold device at Vbias for time tstab
            par.tmesh_type = 1;
            par.t0 = 0;
            par.tmax = tstab;
            par.tpoints = 100;
        
            par.V_fun_type = 'constant';
            par.V_fun_arg(1) = Vbias(i) + 0.01;
            
            disp(['Stabilising solution at ' num2str(Vbias(i) + 0.01) ' V'])

            SaPsol{i,1} = df(sol, par);
    
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
            
                %Hold device at Vbias for time tstab
                par.tmesh_type = 1;
                par.t0 = 0;
                par.tmax = tstab;
                par.tpoints = 100;
            
                par.V_fun_type = 'constant';
                par.V_fun_arg(1) = Vbias(i) - 0.01;
                
                disp(['Stabilising solution at ' num2str(Vbias(i) - 0.01) ' V'])

                SaPsol{i,1} = df(sol, par);
        
                Vbias(i) = Vbias(i) - 0.01;
            catch
            warning(['Could not stabilise device at ' num2str(Vbias(i)) ' V'])
            Vbias(i) = 100;
            end
        end 
    end
end    

%% Perform the Pulsed JV for each Vbias 
for i = 1:num_bias
    if Vbias(i) ~= 100

        disp(['Starting SaP for Vstab = ' num2str(Vbias(i)) ' V'])    

        for j = 1:num_pulses
            par = SaPsol{i,1}.par;
        
            %turn off ion motion for the duration of the pulse
            %assuming that this is a valid assumption
            par.mobseti = 0;
        
            %ramp voltage up to the applied voltage over time tramp 
            %NB, this is t_ramp, not tramp...
            par.tmesh_type = 1;
            par.t0 = 0;
            par.tmax = tramp;
            par.tpoints = 100;
        
            par.V_fun_type = 'sweep';
            par.V_fun_arg(1) = Vbias(i);
            par.V_fun_arg(2) = Vpulse(j);
            par.V_fun_arg(3) = tramp;

            try
                sol = df(SaPsol{i,1}, par);
            catch
                warning(['Could not ramp to voltage for Vpulse = ' num2str(Vpulse(j)) ' V'])
                sol = 0;
            end
            
            if not(isa(sol, 'struct'))
                SaPsol{i,j+1}.Jpulse = 0;
            elseif isa(sol, 'struct')
        
                par = sol.par;
            
                %turn off ion motion for the duration of the pulse
                %assuming that this is a valid assumption
                par.mobseti = 0;
            
                %Hold device at Vpulse for tsample 
                %On paper they say this is 1e-3 seconds after the pulse applied 
                par.tmesh_type = 1;
                par.t0 = 0;
                par.tmax = tsample;
                par.tpoints = 100;
            
                par.V_fun_type = 'constant';
                par.V_fun_arg(1) = Vpulse(j);
                
                disp(['Vpulse = ' num2str(Vpulse(j)) ' V'])

                try
                    SaPsol{i,j+1} = df(sol, par);
                    SaPsol{i,j+1}.Jpulse = dfana.calcJ(SaPsol{i,j+1}).tot(end,1);
                catch
                    warning(['Could not solve Vpulse = ' num2str(Vpulse(j)) ' V'])
                    SaPsol{i,j+1}.Jpulse = 0;
                end 
            end
            
        end 
    else
        for j = 1:num_pulses
            SaPsol{i,j}.Jpulse = 0;
        end
    end
end