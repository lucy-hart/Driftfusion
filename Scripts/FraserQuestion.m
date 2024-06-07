par=pc('Input_files/FraserTest.csv');
% par=pc('Input_files/TiO2_MAPI_Spiro_TestSaP_3.csv');
eqm_QJV = equilibrate(par);

JV_sol_ion = doCV(eqm_QJV.ion, 1, -0.2, 1.3, -0.2, 1e-4, 1, 301);
Plot_Current_Contributions(JV_sol_ion) 
stats_ion = CVstats(JV_sol_ion)