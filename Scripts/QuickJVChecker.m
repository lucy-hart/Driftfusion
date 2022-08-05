%par = pc('Input_files/PTAA_MAPI_NegOffset.csv');
par = pc('Input_files/PTAA_MAPI_NoOffset.csv');
%par = pc('Input_files/PTAA_MAPI_PosOffset.csv');

par.RelTol_vsr = 0.1;
par = refresh_device(par);

eqm_QJV = equilibrate(par);

%%
CV_sol_ion = doCV(eqm_QJV.ion, 1.1, -0.3, 1.2, -0.3, 10e-3, 1, 301);
CV_sol_el = doCV(eqm_QJV.el, 1.1, -0.3, 1.2, -0.3, 10e-3, 1, 301);

%%
Plot_Current_Contributions(CV_sol_ion)
stats_ion = CVstats(CV_sol_ion)
stats_el = CVstats(CV_sol_el)


