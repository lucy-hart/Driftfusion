par = pc('Input_files/3_layer_test.csv');
eqm = equilibrate(par);
CV_solution = doCV(eqm.ion, 1, 0, 1.0, 0, 10e-3, 1, 211);
dfplot.JtotVapp(CV_solution,0)
stats = CVstats(CV_solution);
Plot_Current_Contributions(CV_solution)