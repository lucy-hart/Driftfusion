%par = pc('Input_files/PTAA_MAPI_Kloc6_v2.csv');
%par = pc('Input_files/PTAA_MAPI_PCBM_v2.csv');
par = pc('Input_files/PTAA_MAPI_ICBA_v2.csv');

eqm = equilibrate(par);
CV_sol_ion = doCV(eqm.ion, 1.2, -0.2, 1.2, -0.2, 1e-3, 1, 281);
CV_sol_dark = doCV(eqm.ion, 0, -0.2, 1.2, -0.2, 1e-3, 1, 281);
Plot_Current_Contributions(CV_sol_ion)
stats = CVstats(CV_sol_ion)

e = 1.6e-19;
x = CV_sol_ion.par.x_sub;
gxt = dfana.calcg(CV_sol_ion);
J_gen = e*trapz(x, gxt(1,:))';

%% Calculate series resistance (using sigma = n*e*mu)
% Can just use electron population and mobility as charge conservation
% means that mu_n*n = mu_p*p
% This is spurious but I am proud of it and may suggest that
% Phil add it to dfana

rho_series_el = ((e*CV_sol_ion.par.dev.mu_n.*CV_sol_ion.u(:,:,2))'.^-1)/100; %Convert cm to m
r_series_el = trapz(CV_sol_ion.x, rho_series_el, 1)'/1e6; %x given in um

%% Add a shunt resistance (lab book 21/12/21)
%Kloc6 = 2600 Ohms, PCBM = 35000 Ohms, ICBA = 4700 Ohms
Area = 0.045e-4;
R_shunt = 4700;
r_shunt = Area*R_shunt;
J_PV = dfana.calcJ(CV_sol_ion).tot;
Vapp = dfana.calcVapp(CV_sol_ion)';

J_out = ((Vapp./(r_shunt))*1e-4) + J_PV(:,1); %convert to A/cm^-2

%Calculate JV stats (use forwards scan)
Jsc_f = interp1(Vapp(1:141), J_out(1:141), 0, 'linear');
Voc_f = interp1(J_out(1:141), Vapp(1:141), 0, 'linear');
Pin = dfana.calcPin(CV_sol_ion);
pow_f = J_out(1:141).*Vapp(1:141);
mpp_f = abs(min(pow_f));
efficiency_f = 100*(mpp_f/Pin);
mppV_f = Vapp(-pow_f == mpp_f);
FF_f = -mpp_f/(Jsc_f*Voc_f);

%Plotting for a sanity check
figure(666)
plot(Vapp, J_PV, 'red')
hold on
plot(Vapp, J_out, 'blue')
hold off
xlim([0, 1.3])
ylim([-0.03, 0])

%% Other things I tried

%Differential resistance
%J_dark = dfana.calcJ(CV_sol_dark).tot;
%dV = Vapp(2)-Vapp(1);
%r_grad_dark = (abs(gradient(J_dark(:,1)*1e4, dV))).^-1;
%r_grad_light = (abs(gradient((J_PV(:,1)+stats.Jsc_f)*1e4, dV))).^-1;

%A circuit model (definitely didn't work, I somehow improved Voc...)
%r_diode = abs((Vapp./J_dark(:,1)))*1e-4;
%r_PV = r_diode + r_series_el;

%r_T = (r_shunt*r_PV)./(r_shunt+r_PV);
     