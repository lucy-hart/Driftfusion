par = pc('Input_files/PTAA_MAPI_NoOffset.csv');
eqm_QJV = equilibrate(par);
CV_sol_ion = doCV(eqm_QJV.ion, 1.1, -0.2, 1.2, -0.2, 1e-4, 1, 281);

%%
par_2traps = pc('Input_files/PTAA_MAPI_NoOffset_twotraptest.csv');
eqm_QJV_2traps = equilibrate(par_2traps);
CV_sol_ion_2traps = doCV(eqm_QJV_2traps.ion, 1.1, -0.2, 1.2, -0.2, 1e-4, 1, 281);

%%
Plot_Current_Contributions(CV_sol_ion) 
Plot_Current_Contributions(CV_sol_ion_2traps)
stats_ion = CVstats(CV_sol_ion);
stats_ion_2traps = CVstats(CV_sol_ion_2traps);



