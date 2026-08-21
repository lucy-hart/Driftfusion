par = pc('Input_files/LinneaDonorPolymers.csv');
eqm = equilibrate(par,1);

JVsol = doCV(eqm.el, 0, -5, 1.0, -5, 1, 0.5, 101);

Vapp = dfana.calcVapp(JVsol);
J = dfana.calcJ(JVsol);

figure(666)
plot(Vapp, J.tot(:,1))