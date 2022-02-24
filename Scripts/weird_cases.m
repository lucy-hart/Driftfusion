Ncat = 1;
y_var = 1;

J_ion = dfana.calcJ(solCV_ion{y_var, Ncat});
J_el = dfana.calcJ(solCV_el{y_var, Ncat});
Vapp = dfana.calcVapp(solCV_ion{y_var, Ncat});

figure(008)
plot(Vapp, J_el.tot(:,1), 'blue', Vapp, J_ion.tot(:,1), 'red')

%%
Ncat = 5;
y_var = 4;

solCV_ion_single = doCV(soleq{y_var, Ncat}.ion, 1, -0.1, 1.2, -0.1, 1e-3, 1, 241);
solCV_el_single = doCV(soleq{y_var, Ncat}.el, 1, -0.1, 1.2, -0.1, 1e-4, 1, 241);

