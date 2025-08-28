par = pc('Input_files/SAM_MAFACsPbIBr_C60.csv');
par = refresh_device(par);
eqm = equilibrate(par);
%sol_light = lightonRs(eqm.ion,1,1,1,0,100);
JVsol = doCV(eqm.ion, 1, -0.2, 1.2, -0.2, 10e-3, 1, 281);
Plot_Current_Contributions(JVsol)