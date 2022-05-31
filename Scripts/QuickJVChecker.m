par = pc('Input_files/PTAA_MAPI_Kloc6_v3.csv');
%par = pc('Input_files/PTAA_MAPI_PCBM_v3.csv');
%par = pc('Input_files/PTAA_MAPI_IPH_v3.csv');
%par = pc('Input_files/PTAA_MAPI_ICBA_v3.csv');
par.RelTol_vsr = 0.1;
par = refresh_device(par);

eqm_QJV = equilibrate(par);
CV_sol_ion = doCV(eqm_QJV.ion, 1.1, -0.3, 1.2, -0.3, 1e-3, 1, 321);
%CV_sol_el = doCV(eqm_QJV.el, 1.1, -0.3, 1.3, -0.3, 1e-3, 1, 321);
Plot_Current_Contributions(CV_sol_ion)
stats = CVstats(CV_sol_ion)

%%
num_start = sum(CV_sol_ion.par.layer_points(1:2))+1;
num_stop = num_start + CV_sol_ion.par.layer_points(3)-1;
x = CV_sol_ion.par.x_sub;
d = CV_sol_ion.par.d(3);
[~, ~, Efn, Efp] = dfana.calcEnergies(CV_sol_ion);
QFLS_SC = trapz(x(num_start:num_stop), Efn(31, num_start:num_stop)-Efp(31,num_start:num_stop))/d;
  
%% Calculate QFLS at OC point
%Need to find time point solutuion is evaluated which is closest to Voc first
Vapp = dfana.calcVapp(CV_sol_ion);
Voc = stats.Voc_f;
OC_time = find(abs(Voc-Vapp) == min(abs(Voc-Vapp)),1);
    
[~, ~, Efn, Efp] = dfana.calcEnergies(CV_sol_ion);
QFLS_OC = trapz(x(num_start:num_stop), Efn(OC_time, num_start:num_stop)-Efp(OC_time,num_start:num_stop))/d;
  
%% Find 'Figure of Merit'

Delta_mu = (QFLS_OC-QFLS_SC)*1000;
QFLS_Loss = (QFLS_OC-Voc)*1000;


