%par = pc('Input_files/PTAA_MAPI_Kloc6_v2.csv');
par = pc('Input_files/PTAA_MAPI_PCBM_v2.csv');
%par = pc('Input_files/PTAA_MAPI_ICBA_v2.csv');

eqm_QJV = equilibrate(par);
CV_sol_ion = doCV(eqm_QJV.ion, 1.15, -0.3, 1.3, -0.3, 1e-3, 1, 321);
CV_sol_el = doCV(eqm_QJV.el, 1.15, -0.3, 1.3, -0.3, 1e-3, 1, 321);
Plot_Current_Contributions(CV_sol_ion)
stats = CVstats(CV_sol_ion)

