

%par = pc('Input_files/pedotpss_mapi_kloc6.csv');
par = pc('Input_files/pedotpss_mapi_icba.csv');
%par = pc('Input_files/pedotpss_mapi_pcbm.csv');

eqm = equilibrate(par);
CV_sol_ion = doCV(eqm.ion, 1, -0.3, 1.3, -0.3, 10e-3, 1, 241);
dfplot.JtotVapp(CV_sol_ion,0)