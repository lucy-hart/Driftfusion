function SaPsol = doSaP(sol_ini, Vbias, Vpulse, tpulse, tcycle, tstab, light_intensity)
% Performs a simulation of a stabilise and pulse (SaP) measurement
% Input arguments:
% SOL_INI = solution containing intitial conditions (dark eqm device)
% VBIAS = array of stabilisation biases to sample at
% VPULSE = array of voltages to sample at in the pulsed JV
% TPULSE = the length of the voltage pulse
% TCYCLE = the length of the cycle
% TSTAB = time for which the device is stabilised at Vbias
% LIGHT_INTENSITY = light intensity at which the measurement is done
% 
%

%% Start Code
disp('Starting SaP')

num_bias = length(Vbias);
num_pulses = length(Vpulse);

duty_cycle = 100*(tpulse/tcycle);

%Each row in the solution structure is a different prebais
%The first column is the solution after being held at the relevant Vbias
%for time tstab
%The other columns are the resuls of the pulsed JVs
SaPsol = cell(num_bias, num_pulses+1);

%% Go to correct light intensity
%Funtion assumes you are using the first light source currently 
sol_ill = changeLight(sol_ini, light_intensity, 0, 1);

%% Genrate solutions after voltage stabilisation period
for i = 1:num_bias
    
    par = sol_ill.par;

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

end    

%% Perform the Pulsed JV for each Vbias 
for i = 1:num_bias
    disp(['Starting SaP for Vstab = ' num2str(Vbias(i)) ' V'])
    par = SaPsol{i,1}.par;

    %turn off ion motion for the duration of the pulse
    %assuming that this is a valid assumption
    par.mobseti = 0;

    %Settings for the voltage Pulse 
    par.tmesh_type = 1;
    par.t0 = 0;
    par.tmax = tcycle;
    par.tpoints = 1e3;

    %smoothed-sqaure defined so that pulse starts at time 0.05*tcycle
    par.V_fun_type = 'smoothed_square';
    par.V_fun_arg(1) = Vbias(i);
    par.V_fun_arg(3) = tcycle;
    par.V_fun_arg(4) = duty_cycle;

    for j = 1:num_pulses

        par.V_fun_arg(2) = Vpulse(j);

        disp(['Vpulse = ' num2str(Vpulse(j)) ' V'])
        SaPsol{i,j+1} = df(SaPsol{i,1}, par);

    end
end