par = pc('Input_files/3_layer_test.csv');
eqm = equilibrate(par);
CV_solution = doCV(eqm.ion, 1, 0, 1.1, 0, 10e-3, 1, 211);

%% Plots
dfplot.JtotVapp(CV_solution,0)
Plot_Current_Contributions(CV_solution)
ylim([-25e-3, 25e-3])
stats = CVstats(CV_solution);