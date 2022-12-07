%Code to simulate TPV/TPC measurements

%Excitation provided by a pulsed Continuum Minilite Nd:YAG laser at 532 nm with a pulse width of <10 ns 
%Use time base of 1e-5s to start with 

%Should check what dwell time is after changing to a new illumination
%For now will just let device equilibriate between measurements 

%% Read in data files 
par_pcbm = pc('Input_files/PTAA_MAPI_PCBM_v5.csv');
par_pcbm_LowerMob = pc('Input_files/PTAA_MAPI_PCBM_LowerMob.csv');
par_pcbm_LowestMob = pc('Input_files/PTAA_MAPI_PCBM_LowestMob.csv');
par_pcbm_HigherLUMO = pc('Input_files/PTAA_MAPI_PCBM_HigherLUMO.csv');
par_pcbm_LowerLUMO = pc('Input_files/PTAA_MAPI_PCBM_LowerLUMO.csv');

devices = {par_pcbm, par_pcbm_LowerMob, par_pcbm_LowestMob, par_pcbm_LowerLUMO, par_pcbm_HigherLUMO};
num_devices = length(devices);
%num_devices = 1;
SArea = 0.045;

%Set laser parameters
for i = 1:num_devices
    par = devices{i};
    par.light_source2 = 'laser';
    par.laser_lambda2 = 532;
    par.pulsepow = 1;
    par.RelTol_vsr = 0.1;
    devices{i} = refresh_device(par);
end

%% Find eqm solutions 
eqm_solutions_dark = cell(1,num_devices);
for i = 1:num_devices
    eqm_solutions_dark{i} = equilibrate(devices{i});
end

%% Do TPC
TPC_solutions_illuminated = cell(num_devices);
TPC_solutions_dark = cell(num_devices);
J_transient_ill = cell(num_devices);
J_transient_dark = cell(num_devices);
DeltaQ = zeros(num_devices);
%Add 50 ohm resisistor for TPC transient
%Rs is actually R*device area - based on looking at how it comes into
%setting the BCs in df (Vres = J*Rs)
%Device area in cm2 
pulse_int = 10;

%Got numerical ringing when I tried to use doLightPulse for the entire TPC
%transient. Split solutuion into two parts, one with the laser pulse and
%one after the pulse has finished. 
for i = 1:num_devices
    eqm_solutions_dark{i}.ion.par.Rs = 50*SArea;
    FiftyOhmSolution = stabilize(eqm_solutions_dark{i}.ion);
    TPC_solutions_illuminated{i} = doLightPulse(FiftyOhmSolution, pulse_int, 10e-9, 1000, 100, 1, 1);
    TPC_solutions_dark{i} = changeLight(TPC_solutions_illuminated{i}, 0, 5e-6, 1);
    J_transient_ill{i} = -dfana.calcJ(TPC_solutions_illuminated{i}).tot;
    J_transient_dark{i} = -dfana.calcJ(TPC_solutions_dark{i}).tot;
    DeltaQ(i) = trapz(TPC_solutions_illuminated{i}.t, J_transient_ill{i}(:,1),1);
end

figure('Name', 'TPC')
for i = 1:num_devices
    hold on 
    plot(TPC_solutions_illuminated{i}.t, J_transient_ill{i}(:,end))
    plot(TPC_solutions_dark{i}.t + 10e-9, J_transient_dark{i}(:,end))
end
hold off 
xlabel('Time (s)')
ylabel('Current Density (mAcm^{-2})')
%% Do TPV
%Make light source 1 AM1.5G for the bias light
for i = 1:num_devices
    par = devices{i};
    par.light_source1 = 'AM';
    devices{i} = refresh_device(par);
end

num_samples = 10;
TPV_solutions = cell(num_devices, num_samples);
DeltaV = cell(num_devices, num_samples);
intensities = logspace(log10(0.5), log10(5), num_samples);
%Duty is the percentage of t max than the laser in on for.
tmax = 7e-6;
tlaser = 10e-9;
tduty = 100*(tlaser/tmax);
pulse_int = 1;
for i = 1:num_devices
    for j = 1:num_samples
        TPV_solutions{i,j} = doTPV(eqm_solutions_dark{i}.ion, intensities(j), 2, 1, 1e6, pulse_int, tmax, 1000, tduty);
        DeltaV{i,j} = dfana.deltaVt(TPV_solutions{i,j}, 1, end);
    end
end

figure('Name', 'TPV')
num = 1;
for j = 1:num_samples
    hold on 
    plot(TPV_solutions{num,j}.t, Delta_V{num,j})
end
hold off 
xlabel('Time (s)')
ylabel('Transient Photovoltage (V)')
