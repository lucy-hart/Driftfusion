par1=pc('Input_files/FaCs-PCBM-Charlie.csv');
par2=pc('Input_files/FaCs-Charlie.csv');
% par1.z_c = -1;
% par2.z_c = -1;
eqm_calib = equilibrate(par1);


%% Do calicration for voltages

JV_w_PCBM = doCV(eqm_calib.ion, 0, -0.2, 1.2, -0.2, 1e-4, 0.5, 141);
p_sample = sum(JV_w_PCBM.par.layer_points(1:2));
V_pero = JV_w_PCBM.u(:,p_sample,1);
V_app = dfana.calcVapp(JV_w_PCBM);
Vbi_eff = V_pero(21);

par2.Phi_right = par2.Phi_left + Vbi_eff;
eqm_JV = equilibrate(par2);
%%
suns = 1;
V_bias = Vbi_eff + 0.2;
V_max = V_bias;
V_min = -0.1;
scan_rate = 10e-3;
deltaV = V_max - V_min;
tmax = deltaV/scan_rate;
biased_eqm_ion = genVappStructs(eqm_JV.ion, V_bias, 0);
illuminated_sol_ion = changeLight(biased_eqm_ion, suns, 0, 1);
JV_sol_ion_rev = VappFunction(illuminated_sol_ion, 'sweep', [V_max, V_min, tmax], tmax, 200*(V_max-V_min)+1, 0);
JV_sol_ion_fw = VappFunction(JV_sol_ion_rev, 'sweep', [V_min, V_max, tmax], tmax, 200*(V_max-V_min)+1, 0);

%% Plot JVs
v_fw = dfana.calcVapp(JV_sol_ion_fw);
v_pero_fw = JV_sol_ion_fw.u(:,end,1);
v_fw_eff = interp1(V_pero, V_app', v_pero_fw);

v_rev = dfana.calcVapp(JV_sol_ion_rev);
v_pero_rev = JV_sol_ion_rev.u(:,end,1);
v_rev_eff = interp1(V_pero, V_app', v_pero_rev);

figure('Name', 'JVPlot', 'Position', [100 100 1250 1250])
hold on
xline(0, 'black', 'LineWidth', 1)
yline(0, 'black', 'LineWidth', 1)
j_fw = dfana.calcJ(JV_sol_ion_fw).tot(:,1);
j_rev = dfana.calcJ(JV_sol_ion_rev).tot(:,1);
plot(v_fw(1:end), j_fw(1:end)*1000, 'color', [0.4660 0.6740 0.1880], 'LineWidth', 3) 
plot(v_rev(1:end), j_rev(1:end)*1000, 'color', [0.4660 0.6740 0.1880], 'LineWidth', 3, 'LineStyle', '--') 
plot(v_fw_eff(1:end), j_fw(1:end)*1000, 'color', [0 0.4470 0.7410], 'LineWidth', 3) 
plot(v_rev_eff(1:end), j_rev(1:end)*1000, 'color', [0 0.4470 0.7410], 'LineWidth', 3, 'LineStyle', '--') 
hold off

box on 
set(gca, 'FontSize', 25)
xlim([-0.15, 1.2])
ylim([-25,5])
%legend({'','','Mobile Ions','No Mobile Ions'}, 'Location', 'northwest', 'FontSize', 30)
xlabel('Voltage(V)', 'FontSize', 30)
ylabel('Current Density (mAcm^{-2})', 'FontSize', 30)
ax1 = gcf;
