% parC60 = pc('Input_files/HTL_MAPI_C60_DavideValues.csv');
parC60 = pc('Input_files/HTL_MAPI_C60_DavideValues_dark.csv');
parPM6 = pc('Input_files/HTL_MAPI_PM6Y6_DavideValues.csv');
% parPM7 = pc('Input_files/HTL_MAPI_PM7Y6_DavideValues.csv');
% parPCE12 = pc('Input_files/HTL_MAPI_PCE12Y6_DavideValues.csv');

eqm_QJV_C60 = equilibrate(parC60);
% eqm_QJV_PM6 = equilibrate(parPM6);
%% 
%Check if the dark eqm solution has a non-zero current 
%i.e., calibrate for errors in numerical integration

voltage_ar = linspace(-0.5, 0.5, 21);
Jdark = doDarkJV(eqm_QJV_C60.ion, voltage_ar, 1);

%%
figure('Name', 'Dark JV')
plot(voltage_ar, Jdark.Jvalue)

xlabel('Voltage (V)')
xlim([voltage_ar(1), voltage_ar(end)])
ylabel('Current Density (mA cm^{-2})')
ylim([-1e-10, 1e-10])