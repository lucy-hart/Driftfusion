par_10 = pc('Input_files/ptaa_mapi_pcbm_with_vsr.csv');
eqm = equilibrate(par_10);
CV_solution = doCV(eqm.ion, 1, 0, 1.3, 0, 10e-3, 1, 241);
Plot_Current_Contributions(CV_solution)