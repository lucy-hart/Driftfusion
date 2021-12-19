%par = pc('Input_files/PTAA_MAPI_Kloc6_v2.csv');
%par = pc('Input_files/PTAA_MAPI_PCBM_v2.csv');
par = pc('Input_files/PTAA_MAPI_ICBA_v2.csv');

eqm = equilibrate(par);
CV_sol_ion = doCV(eqm.ion, 1.15, -0.3, 1.3, -0.3, 1e-3, 1, 341);
dfplot.JtotVapp(CV_sol_ion,0)
Plot_Current_Contributions(CV_sol_ion)
legend('data', 'Location', 'northwest')
CVstats(CV_sol_ion)

%% Add shunt resistance 

J_tot = dfana.calcJ(CV_sol_ion).tot;
V = dfana.calcVapp(CV_sol_ion)';
R_shunt = 20;
J_ext = J_tot(:,1) + V/R_shunt;
figure(7)
plot(V, J_ext)
hold on
plot(V, J_tot(:,1))
hold off

     