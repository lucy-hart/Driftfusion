%par=pc('Input_files/SAM_MAFACsPbIBr_C60.csv');
par=pc('Input_files/SAMFree_MAFACsPbIBr_C60.csv');
% par.Phi_left = -4.7;
% par.sn_l = 100;
eqm_QJV = equilibrate(par);
sol_bias = genVappStructs(eqm_QJV.ion, 0.8, 1);
sol_bias.par.mobseti = 0;

% JV_sol_el = doCV(eqm_QJV.el, suns, -0.2, 1.25, -0.2, 1, 1, 291);
JV_sol_el = doCV(sol_bias, suns, -0.2, 1.25, -0.2, 1, 1, 291);
JV_sol_ion = doCV(eqm_QJV.ion, suns, -0.2, 1.25, -0.2, 1e-3, 1, 291);

Plot_Current_Contributions(JV_sol_ion)
Plot_Current_Contributions(JV_sol_el)

%% Plot JVs
figure('Name', 'JVPlot', 'Position', [100 100 1250 1250])
colors_JV = {[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],[0 0.4470 0.7410],[0.9290 0.6940 0.1250]};

v_el = dfana.calcVapp(JV_sol_el)';
v_ion = dfana.calcVapp(JV_sol_ion)';

hold on
xline(0, 'black', 'LineWidth', 1)
yline(0, 'black', 'LineWidth', 1)
j_el = dfana.calcJ(JV_sol_el).tot(:,1);
j_ion = dfana.calcJ(JV_sol_ion).tot(:,1);  
plot(v_ion(1:end), j_ion(1:end)*1000, 'color', [0.4660 0.6740 0.1880], 'LineWidth', 3) 
plot(v_el(1:end), j_el(1:end)*1000, 'color', [0 0.4470 0.7410], 'LineWidth', 3)
hold off

box on 
set(gca, 'FontSize', 25)
xlim([-0.15, 1.2])
ylim([-25,5])
legend({'','','With Mobile Ions','Without Mobile Ions'}, 'Location', 'northwest', 'FontSize', 30)
xlabel('Voltage(V)', 'FontSize', 30)
ylabel('Current Density (mAcm^{-2})', 'FontSize', 30)
ax1 = gcf;