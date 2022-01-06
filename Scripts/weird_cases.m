V_bi = 1;
surf_vel = 1;

J_ion = dfana.calcJ(solCV_ion{surf_vel, V_bi});
J_el = dfana.calcJ(solCV_el{surf_vel, V_bi});
Vapp = dfana.calcVapp(solCV_ion{surf_vel, V_bi});

figure(008)
plot(Vapp, J_el.tot(:,1), 'blue', Vapp, J_ion.tot(:,1), 'red')

%% SS solutions under illumination

el_light = changeLight(soleq{surf_vel, V_bi}.el, 1, 0);
ion_light = changeLight(soleq{surf_vel, V_bi}.ion, 1, 0);



