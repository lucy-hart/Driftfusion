%Check that 3 layer symmetric device actually works like I expect it to. 

par = pc('Input_files/3_layer_test_symmetric.csv');

eqm = equilibrate(par);

solCV_ion = doCV(eqm.ion, 1, -0.1, 1.6, -0.1, 1e-3, 1, 241);
solCV_el = doCV(eqm.el, 1, -0.1, 1.6, -0.1, 1e-3, 1, 241);
v = dfana.calcVapp(solCV_ion);

%% Make Plot
figure(1)
plot(v, dfana.calcJ(solCV_ion).tot(:,1), 'r',  v, dfana.calcJ(solCV_el).tot(:,1), 'b')
xlim([-0.1, 1.7])
ylim([-0.025, 0.01])

%% Get Stats
stats_ion = CVstats(solCV_ion);
stats_el = CVstats(solCV_el);