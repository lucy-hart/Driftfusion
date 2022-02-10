%Calculate QFLS for Weidong ETL stuff. Run this after Weidong_ETL

%% Calculate QFLS at SC - do el and ion cases
num_start = sum(CV_solutions_ion{1}.par.layer_points(1:2))+1;
num_stop = num_start + CV_solutions_ion{1}.par.layer_points(3)-1;
x = CV_solutions_ion{1}.par.x_sub;
d = CV_solutions_ion{1}.par.d(3);
QFLS_SC_ion = zeros(1,3);
QFLS_SC_el = zeros(1,3);
for y=1:3
    [Ecb_ion, Evb_ion, Efn_ion, Efp_ion] = dfana.calcEnergies(CV_solutions_ion{y});
    [Ecb_el, Evb_el, Efn_el, Efp_el] = dfana.calcEnergies(CV_solutions_el{y});
    QFLS_SC_ion(y) = trapz(x(num_start:num_stop), Efn_ion(31, num_start:num_stop)-Efp_ion(31,num_start:num_stop))/d;
    QFLS_SC_el(y) = trapz(x(num_start:num_stop), Efn_el(31,num_start:num_stop)-Efp_el(31,num_start:num_stop))/d;
end
%% Calculate QFLS at OC
%Need to find time point solutuion is evaluated at closest to Voc first
Voc_ion = zeros(1,3);
Voc_el = zeros(1,3);
OC_time_ion = zeros(1,3);
OC_time_el = zeros(1,3);
QFLS_OC_ion = zeros(1,3);
QFLS_OC_el = zeros(1,3);

for z=1:3
    V_temp = dfana.calcVapp(CV_solutions_ion{z});
    Voc_ion(z) = CVstats(CV_solutions_ion{z}).Voc_f;
    OC_time_ion(z) = find(abs(Voc(z)-V_temp) == min(abs(Voc(z)-V_temp)),1);
    Voc_el(z) = CVstats(CV_solutions_el{z}).Voc_f;
    OC_time_el(z) = find(abs(Voc(z)-V_temp) == min(abs(Voc(z)-V_temp)),1);
end

for w=1:3
    [Ecb_ion, Evb_ion, Efn_ion, Efp_ion] = dfana.calcEnergies(CV_solutions_ion{w});
    [Ecb_el, Evb_el, Efn_el, Efp_el] = dfana.calcEnergies(CV_solutions_el{w});
    QFLS_OC_ion(w) = trapz(x(num_start:num_stop), Efn_ion(OC_time_ion(w), num_start:num_stop)-Efp_ion(OC_time_ion(w),num_start:num_stop))/d;
    QFLS_OC_el(w) = trapz(x(num_start:num_stop), Efn_el(OC_time_el(w), num_start:num_stop)-Efp_el(OC_time_el(w),num_start:num_stop))/d;
end

%% Find (and print)'Figure of Merit'

Delta_mu_ion = (QFLS_OC_ion-QFLS_SC_ion)*1000
Delta_mu_el = (QFLS_OC_el-QFLS_SC_el)*1000
QFLS_Loss = (QFLS_OC_ion-Voc_ion)*1000
