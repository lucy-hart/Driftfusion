par = pc('Input_files/PTAA_MAPI_PCBM_v5.csv');
%par = pc('Input_files/PTAA_MAPI_PCBM_HigherLUMO.csv');

par.light_source1 = 'laser';
par.laser_lambda1 = 532;
par.pulsepow = 62;
par.RelTol_vsr = 0.1;
par = refresh_device(par);

eqm_QJV = equilibrate(par);

%%
CV_sol_ion = doCV(eqm{1}.ion, LightInt(5), -0.3, 1.2, -0.3, 10e-3, 1, 301);
%CV_sol_el = doCV(eqm_QJV.el, 1.15, -0.3, 1.2, -0.3, 10e-3, 1, 301);

%%
Plot_Current_Contributions(CV_sol_ion,0)
stats = CVstats(CV_sol_ion)

%%
num_start = sum(CV_sol_ion.par.layer_points(1:2))+1;
num_stop = num_start + CV_sol_ion.par.layer_points(3)-1;
x = CV_sol_ion.par.x_sub;
d = CV_sol_ion.par.d(3);
[~, ~, Efn, Efp] = dfana.calcEnergies(CV_sol_ion);
QFLS_SC = trapz(x(num_start:num_stop), Efn(291, num_start:num_stop)-Efp(291,num_start:num_stop))/d;
  
%% Calculate QFLS at OC point
%Need to find time point solutuion is evaluated which is closest to Voc first
Vapp = dfana.calcVapp(CV_sol_ion);
Voc = stats.Voc_r;
OC_time = find(abs(Voc-Vapp) == min(abs(Voc-Vapp)),1);
    
[~, ~, Efn, Efp] = dfana.calcEnergies(CV_sol_ion);
QFLS_OC = trapz(x(num_start:num_stop), Efn(OC_time, num_start:num_stop)-Efp(OC_time,num_start:num_stop))/d;
  
%% Find 'Figure of Merit'

Delta_mu = (QFLS_OC-QFLS_SC)*1000;
QFLS_Loss = (QFLS_OC-Voc)*1000;


