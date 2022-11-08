par = pc("Input_files/ptaa_mapi_pcbm.csv");
eqm = equilibrate(par);
CV_sol_ion = doCV(eqm.ion, 1, -0.3, 1.3, -0.3, 1e-3, 1, 321);
CV_sol_el = doCV(eqm.el, 1, -0.3, 1.3, -0.3, 1e-3, 1, 321);

Vapp = dfana.calcVapp(CV_sol_ion);
Jion = dfana.calcJ(CV_sol_ion);
Jel = dfana.calcJ(CV_sol_el);

figure(666)
plot(Vapp, Jion.tot(:,1), Vapp, Jel.tot(:,1))