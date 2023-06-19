filenames = ['Input_files/HTL_MAPI_C60_FrancescoValues.csv',...
            'Input_files/HTL_MAPI_PM6Y6_FrancescoValues.csv',...            
            'Input_files/HTL_MAPI_Y6_FrancescoValues.csv'];

num_devices = length(filenames);
devices = cell(1, num_devices);
for i = 1:num_devices
    devices{i} = pc(filenames(i));
end

%%
run_C60 = 1;

if run_C60 == 1
    eqm_QJV_C60 = equilibrate(parC60);
    CV_sol_C60 = doCV(eqm_QJV_C60.ion, 1.0, -0.2, 1.25, -0.2, 25e-3, 1, 291);
    %CV_sol_C60_el = doCV(eqm_QJV_C60.el, 1.0, -0.2, 1.25, -0.2, 25e-3, 1, 291);
    stats_C60 = CVstats(CV_sol_C60);
    %stats_C60_el = CVstats(CV_sol_C60_el);
    J_C60 = dfana.calcJ(CV_sol_C60);
    Plot_Current_Contributions(CV_sol_C60)
end

eqm_QJV_PM6Y6 = equilibrate(parPM6Y6);
CV_sol_PM6Y6 = doCV(eqm_QJV_PM6Y6.ion, 1.0, -0.2, 1.25, -0.2, 25e-3, 1, 291);
%CV_sol_PM6Y6_el = doCV(eqm_QJV_PM6Y6.el, 1.0, -0.2, 1.25, -0.2, 25e-3, 1, 291);
stats_PM6Y6 = CVstats(CV_sol_PM6Y6);
%stats_PM6Y6_el = CVstats(CV_sol_PM6Y6_el);
Plot_Current_Contributions(CV_sol_PM6Y6)

Vapp = dfana.calcVapp(CV_sol_PM6Y6);
J_PM6Y6 = dfana.calcJ(CV_sol_PM6Y6);

%%
if run_C60 == 1
    figure('Name', 'JV Plots', 'Position', [50 50 1000 1000])

    box on
    hold on
    xline(0, 'LineWidth', 2, 'Color', 'black')
    yline(0, 'LineWidth', 2, 'Color', 'black')
    plot(Vapp, 1e3*J_C60.tot(:,1), 'LineWidth', 4, 'Color', [0 0.4470 0.7410])
    plot(Vapp, 1e3*J_PM6Y6.tot(:,1), 'LineWidth', 4, 'Color', [0.4660 0.6740 0.1880])
    hold off

    set(gca, 'FontSize', 25)
    xlabel('Voltage (V)', 'FontSize', 30)
    ylabel('Current Density (mA cm^{-2})', 'FontSize', 30)
    xlim([-0.15, 1.2])
    ylim([-25, 5])
    legend({'', '', 'C_{60}', 'PM6:Y6'}, 'FontSize', 25, 'Location', 'northwest')
end 


