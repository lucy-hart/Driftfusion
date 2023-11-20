par = pc('Input_files/SAM_MAFACsPbIBr_PM6Y6.csv');
% par = pc('Input_files/SAM_MAFACsPbIBr_C60.csv');
eqm = equilibrate(par);

CV_sol_ion = doCV(eqm.ion, 1, -0.3, 1.3, -0.3, 1e-3, 1, 321);
% CV_sol_el = doCV(eqm.el, 1, -0.3, 1.3, -0.3, 1e-3, 1, 321);

Vapp = dfana.calcVapp(CV_sol_ion);
Jion = dfana.calcJ(CV_sol_ion);
% Jel = dfana.calcJ(CV_sol_el);

figure(666)
plot(Vapp, Jion.tot(:,1))
