par=pc('Input_files/1_layer_test.csv');

par.Phi_left = par.EF0-0.1;
par.Phi_right = par.EF0+0.1;
par.sn_l = 1e6;
par.sp_r = 1e6;
par = refresh_device(par);
soleq = equilibrate(par);

solCV_ion = doCV(soleq.ion, 1, -0.2, 1.2, -0.2, 1e-4, 1, 241);
solCV_el = doCV(soleq.el, 1, -0.2, 1.2, -0.2, 1e-4, 1, 241);

J_ion = dfana.calcJ(solCV_ion);
J_el = dfana.calcJ(solCV_el);

figure(008)
plot(dfana.calcVapp(solCV_el), J_el.tot(:,1), 'blue', dfana.calcVapp(solCV_ion), J_ion.tot(:,1), 'red')

el_light = changeLight(soleq.el, 1, 0);
ion_light = changeLight(soleq.ion, 1, 0);


