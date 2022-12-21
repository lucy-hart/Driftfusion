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

LightInt = logspace(log10(0.1), log10(5), num_samples);
TPVsols = cell(num_devices, num_samples);

%% Do TPV Measurement
%Check that the duty cycle is right 

for i = 1:num_devices
    for j = 1:num_samples
    TPVsols{i,j} = doTPV(eqm{i}.ion, 1, -1, 1, 1e6, 5, 9e-6, 1000, 0.1);
    end
end

%Make this so it subtracts the Voc, just want to plot the DeltaV. 
dfplot.Voct(TPVsols{1,1})