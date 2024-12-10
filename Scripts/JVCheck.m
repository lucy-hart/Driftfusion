par = pc('Input_files/EnergyOffsetSweepParameters_v5_undoped_SAM.csv');
eqm = equilibrate(par);
CV_sol_ion = doCV(eqm.ion, 1, -0.3, 1.7, -0.3, 1e-4, 1, 401);
%%
Vapp = dfana.calcVapp(CV_sol_ion);
Jion = dfana.calcJ(CV_sol_ion);
Jel = dfana.calcJ(CV_sol_el);

figure(666)
plot(Vapp, Jion.tot(:,1), Vapp, Jel.tot(:,1))