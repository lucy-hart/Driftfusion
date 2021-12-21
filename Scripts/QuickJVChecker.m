par = pc('Input_files/PTAA_MAPI_Kloc6_v2.csv');
%par = pc('Input_files/PTAA_MAPI_PCBM_v2.csv');
%par = pc('Input_files/PTAA_MAPI_ICBA_v2.csv');

eqm = equilibrate(par);
CV_sol_ion = doCV(eqm.ion, 1.15, -0.3, 1.3, -0.3, 1e-3, 1, 341);
dfplot.JtotVapp(CV_sol_ion,0)
Plot_Current_Contributions(CV_sol_ion)
legend('data', 'Location', 'northwest')
CVstats(CV_sol_ion)

%% Try and add a shunt resistance

Area = 0.045e-4;
R_shunt = 3e6;
J_SS = dfana.calcJ(CV_sol_ion).tot;
Vapp = dfana.calcVapp(CV_sol_ion)';
R_SS = Vapp./(J_SS*Area);
R_T = (R_shunt*R_SS)./(R_shunt+R_SS);
J_out = Vapp./(R_T*Area);
figure(666)
plot(Vapp, J_SS)
hold on
plot(Vapp, J_out)
hold off
legend('Driftfusion', 'With shunt?')





     