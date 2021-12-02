%par = pc('Input_files/PTAA_MAPI_Kloc6_SmallV_BI.csv');
par = pc('Input_files/PTAA_MAPI_PCBM_SmallV_BI.csv');
%par = pc('Input_files/PTAA_MAPI_ICBA_SmallV_BI.csv');

eqm = equilibrate(par);
CV_sol_ion = doCV(eqm.ion, 1, -0.3, 1.3, -0.3, 10e-3, 1, 241);
dfplot.JtotVapp(CV_sol_ion,0)