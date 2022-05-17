%Calculate QFLS for Weidong ETL stuff. Run this after Weidong_ETL

%% Calculate QFLS 
num_start = sum(CV_solutions_ion{1}.par.layer_points(1:2))+1;
num_stop = num_start + CV_solutions_ion{1}.par.layer_points(3)-1;
num_points = length(CV_solutions_ion{1}.t);
x = CV_solutions_ion{1}.par.x_sub;
d = CV_solutions_ion{1}.par.d(3);
QFLS_ion = zeros(num_points,num_devices);
QFLS_el = zeros(num_points,num_devices);
for y=1:num_devices
    [Ecb_ion, Evb_ion, Efn_ion, Efp_ion] = dfana.calcEnergies(CV_solutions_ion{y});
    [Ecb_el, Evb_el, Efn_el, Efp_el] = dfana.calcEnergies(CV_solutions_el{y});
    QFLS_ion(:,y) = trapz(x(num_start:num_stop), Efn_ion(:, num_start:num_stop)-Efp_ion(:,num_start:num_stop),2)/d;
    QFLS_el(:,y) = trapz(x(num_start:num_stop), Efn_el(:,num_start:num_stop)-Efp_el(:,num_start:num_stop),2)/d;
end

QFLS_SC_ion = QFLS_ion(31,:);
QFLS_SC_el = QFLS_el(31,:);

%% Calculate OC time points
%Need to find time point solutuion is evaluated at closest to Voc first
Voc_ion = zeros(1,num_devices);
Voc_el = zeros(1,num_devices);
OC_time_ion = zeros(1,num_devices);
OC_time_el = zeros(1,num_devices);
QFLS_OC_ion = zeros(1,num_devices);
QFLS_OC_el = zeros(1,num_devices);

for z=1:num_devices
    V_temp = dfana.calcVapp(CV_solutions_ion{z});
    Voc_ion(z) = CVstats(CV_solutions_ion{z}).Voc_f;
    OC_time_ion(z) = find(abs(Voc_ion(z)-V_temp) == min(abs(Voc_ion(z)-V_temp)),1);
    QFLS_OC_ion(z) = QFLS_ion(OC_time_ion(z), z);
    Voc_el(z) = CVstats(CV_solutions_el{z}).Voc_f;
    OC_time_el(z) = find(abs(Voc_el(z)-V_temp) == min(abs(Voc_el(z)-V_temp)),1);
    QFLS_OC_el(z) = QFLS_ion(OC_time_el(z), z);
end

%% Find 'Figure of Merit'

Delta_mu_ion = (QFLS_OC_ion-QFLS_SC_ion)*1000;
Delta_mu_el = (QFLS_OC_el-QFLS_SC_el)*1000;
QFLS_Loss = (QFLS_OC_ion-Voc_ion)*1000;

%% Plot QFLS vs Vapp for el vs ion cases

figure(9)
for w=1:num_devices
    plot(V(1:161), QFLS_ion(1:161,w), 'color', colors_JV{w}) 
    hold on 
    semilogy(V(1:161), QFLS_el(1:161,w), '-.', 'color', colors_JV{w})     
    hold on
end
hold off
xlim([0, 1.3])
xlabel('Voltage (V)')
ylabel('QFLS (eV)')
ylim([0.85, 1.2])
if num_devices == 4
    legend({'Kloc-6','', 'PCBM','', 'ICBA','','IPH',''}, 'Location', 'northwest')
elseif num_devices == 3
    legend({'Kloc-6','', 'PCBM','', 'ICBA',''}, 'Location', 'northwest')
end

