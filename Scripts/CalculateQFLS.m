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
OC_time = zeros(1,3);
Voc = zeros(1,3);
QFLS_OC_ion = zeros(1,3);

for z=1:3
    Voc(z) = CVstats(CV_solutions_ion{z}).Voc_f;
    V = dfana.calcVapp(CV_solutions_ion{z});
    OC_time(z) = find(abs(V-Voc(z))== min(abs(V-Voc(z))), 1);
end

for w=1:3
    [Ecb_ion, Evb_ion, Efn_ion, Efp_ion] = dfana.calcEnergies(CV_solutions_ion{w});
    QFLS_OC_ion(w) = trapz(x(num_start:num_stop), Efn_ion(OC_time(w), num_start:num_stop)-Efp_ion(OC_time(w),num_start:num_stop))/d;
end

%% Find (and print)'Figure of Merit'

Delta_mu = (QFLS_OC_ion-QFLS_SC_ion)*1000
QFLS_Loss = (QFLS_OC_ion-Voc)*1000
