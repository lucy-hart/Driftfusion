%% Read in data
par1 = pc('Input_files/HTL_MAPI_NegOffset.csv');
par2 = pc('Input_files/HTL_MAPI_NegOffset.csv');
par3 = pc('Input_files/HTL_MAPI_PosOffset.csv');

pars = {par1, par2, par3};
num_devices = 1;
eqm = cell(1,num_devices);

for i = 1:num_devices
    eqm{i} = equilibrate(pars{i});
end
%% Set up parameters 
num_samples = 10;
num_tpoints = 1000;

LightInt = logspace(log10(0.1), log10(5), num_samples);

TPVsols = cell(num_devices, num_samples);
DeltaVt = zeros(num_devices, num_tpoints, num_samples);
Veqm = zeros(num_devices, num_samples);

TPVsols_noionmob = cell(num_devices, num_samples);
DeltaVt_noionmob = zeros(num_devices, num_tpoints, num_samples);
Veqm_noionmob = zeros(num_devices, num_samples);

%% Do TPV Measurement
%Check that the duty cycle is right 
%Currenly will only work if you run LI_n_Dependence first

for i = 1:num_devices
    for j = 1:num_samples
    TPVsols{i,j} = doTPV(eqm{i}.ion, LightInt(j), -1, 1, 1e6, 1, 9e-6, num_tpoints, 0.1);
    TPVsols_noionmob{i,j} = doTPV(eqm{i}.ion, LightInt(j), -1, 0, 1e6, 1, 9e-6, num_tpoints, 0.1);
    end
end

%% Get Delta V on its own 
for i = 1:num_devices
    for j = 1:num_samples
    Vt = dfana.calcDeltaQFL(TPVsols{i,j});
    DeltaVt(i,:,j) = Vt - Vt(end);
    Veqm(i,j) = Vt(end);
    
    Vt_noionmob = dfana.calcDeltaQFL(TPVsols_noionmob{i,j});
    DeltaVt_noionmob(i,:,j) = Vt_noionmob - Vt_noionmob(end);
    Veqm_noionmob(i,j) = Vt_noionmob(end);
    end
end

%% Plot Delta V vs t
%Make this so it subtracts the Voc, just want to plot the DeltaV. 
figure('Name', 'DeltaV vs t')
t = TPVsols{1,1}.t;

box on
for i = 1:num_devices
    for j = 1:num_samples
        plot(t,DeltaVt(i,:,j))
        hold on
    end
end

%% Plot Delta V vs Voc 
%Make this so it subtracts the Voc, just want to plot the DeltaV. 
figure('Name', 'DeltaV vs Voc')
t = TPVsols{1,1}.t;

box on
for i = 1:num_devices
    for j = 1:num_samples
        semilogy(Veqm(i,j),1/max(DeltaVt(i,:,j)), 'Marker', '+', 'MarkerSize', 10, 'color', 'blue')
        hold on
        semilogy(Veqm_noionmob(i,j),1/max(DeltaVt_noionmob(i,:,j)), 'Marker', '+', 'MarkerSize', 10, 'color', 'red')
    end
end