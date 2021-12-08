par = pc('Input_files/PTAA_MAPI_Kloc6.csv');
%par = pc('Input_files/PTAA_MAPI_PCBM.csv');
%par = pc('Input_files/PTAA_MAPI_ICBA.csv');

eqm = equilibrate(par);
CV_sol_ion = doCV(eqm.ion, 1, -0.3, 1.4, -0.3, 1e-3, 1, 241);
dfplot.JtotVapp(CV_sol_ion,0)
legend('data', 'Location', 'northwest')
