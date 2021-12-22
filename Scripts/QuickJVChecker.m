par = pc('Input_files/PTAA_MAPI_Kloc6_v2.csv');
%par = pc('Input_files/PTAA_MAPI_PCBM_v2.csv');
%par = pc('Input_files/PTAA_MAPI_ICBA_v2.csv');

eqm = equilibrate(par);
CV_sol_ion = doCV(eqm.ion, 1.15, -0.3, 1.3, -0.3, 1e-3, 1, 341);
dfplot.JtotVapp(CV_sol_ion,0)
Plot_Current_Contributions(CV_sol_ion)
legend('data', 'Location', 'northwest')

%% Add a shunt resistance (lab book 21/12/21)

%Calculate the output current
Area = 0.045e-4;
R_shunt = 3e7;
r_shunt = Area*R_shunt;
J_PV = dfana.calcJ(CV_sol_ion).tot;
Vapp = dfana.calcVapp(CV_sol_ion)';
r_SS = Vapp./J_PV(:,1);
r_T = (r_shunt*r_SS)./(r_shunt+r_SS);
J_out = Vapp./(r_T);

%Calculate JV stats (use forwards scan)
Jsc_f = interp1(Vapp(1:171), J_out(1:171), 0, 'linear')
Voc_f = interp1(J_out(1:171), Vapp(1:171), 0, 'linear')
Pin = dfana.calcPin(CV_sol_ion);
pow_f = J_out(1:171).*Vapp(1:171);
mpp_f = abs(min(pow_f));
efficiency_f = 100*(mpp_f/Pin);
mppV_f = Vapp(-pow_f == mpp_f);
FF_f = -mpp_f/(Jsc_f*Voc_f)

%Plotting for a sanity check
figure(666)
plot(Vapp, J_PV, 'red')
hold on
plot(Vapp, J_out, 'blue')
hold off
xlim([0, 1.3])
ylim([-0.03, 0.01])





     