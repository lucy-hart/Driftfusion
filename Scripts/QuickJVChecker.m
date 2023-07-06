parC60 = pc('Input_files/HTL_MAPI_C60_DavideValues.csv');
parPM6 = pc('Input_files/HTL_MAPI_PM6Y6_DavideValues.csv');
%parPM7 = pc('Input_files/HTL_MAPI_PM7Y6_DavideValues.csv');
%parPCE12 = pc('Input_files/HTL_MAPI_PCE12Y6_DavideValues.csv');


%%
run_C60 = 1;
light = 1;
noions = 1;

if light == 1
    suns = 1.0;
    Vmin = -0.2;
elseif light == 0
    suns = 0;
    Vmin = -0.5;
end

if run_C60 == 1
    eqm_QJV_C60 = equilibrate(parC60);
    CV_sol_C60 = doCV(eqm_QJV_C60.ion, suns, Vmin, 1.25, Vmin, 25e-3, 1, 291);
    J_C60 = dfana.calcJ(CV_sol_C60);
    if light == 1
        stats_C60 = CVstats(CV_sol_C60);
        Plot_Current_Contributions(CV_sol_C60)
    end
end

eqm_QJV_PM6 = equilibrate(parPM6);
CV_sol_PM6 = doCV(eqm_QJV_PM6.ion, suns, Vmin, 1.25, Vmin, 25e-3, 1, 291);

%eqm_QJV_PM7 = equilibrate(parPM7);
%CV_sol_PM7 = doCV(eqm_QJV_PM7.ion, suns, Vmin, 1.25, Vmin, 25e-3, 1, 291);

%eqm_QJV_PCE12 = equilibrate(parPCE12);
%CV_sol_PCE12 = doCV(eqm_QJV_PCE12.ion, suns, Vmin, 1.25, Vmin, 25e-3, 1, 291);

if noions == 1
    if run_C60 == 1
        CV_sol_C60_noions = doCV(eqm_QJV_C60.el, suns, Vmin, 1.25, Vmin, 25e-3, 1, 291);
        J_C60_noions = dfana.calcJ(CV_sol_C60_noions);
    end
    CV_sol_PM6_noions = doCV(eqm_QJV_PM6.el, suns, Vmin, 1.25, Vmin, 25e-3, 1, 291);
    J_PM6Y6_noions = dfana.calcJ(CV_sol_PM6_noions);
end

Vapp = dfana.calcVapp(CV_sol_PM6);
J_PM6Y6 = dfana.calcJ(CV_sol_PM6);
%J_PM7Y6 = dfana.calcJ(CV_sol_PM7);
%J_PCE12Y6 = dfana.calcJ(CV_sol_PCE12);

%%
if run_C60 == 1 && light == 1
    figure('Name', 'JV Plots', 'Position', [50 50 1000 1000])

    box on
    hold on
    xline(0, 'LineWidth', 2, 'Color', 'black')
    yline(0, 'LineWidth', 2, 'Color', 'black')
    plot(Vapp, 1e3*J_C60.tot(:,1), 'LineWidth', 4, 'Color', [0 0.4470 0.7410])
    plot(Vapp, 1e3*J_PM6Y6.tot(:,1), 'LineWidth', 4, 'Color', [0.4660 0.6740 0.1880])
    if noions == 1
        plot(Vapp(1:145), 1e3*J_C60_noions.tot(1:145,1), 'LineWidth', 4, 'Color', [0 0.4470 0.7410], 'LineStyle', ':')
        plot(Vapp(1:145), 1e3*J_PM6Y6_noions.tot(1:145,1), 'LineWidth', 4, 'Color', [0.4660 0.6740 0.1880], 'LineStyle', ':')
    end
    %plot(Vapp, 1e3*J_PM7Y6.tot(:,1), 'LineWidth', 4, 'Color', 'red')
    %plot(Vapp, 1e3*J_PCE12Y6.tot(:,1), 'LineWidth', 4, 'Color', 'black')
    hold off

    set(gca, 'FontSize', 25)
    xlabel('Voltage (V)', 'FontSize', 30)
    ylabel('Current Density (mA cm^{-2})', 'FontSize', 30)
    xlim([-0.15, 1.2])
    ylim([-25, 5])
    legend({'', '', ' C_{60}', ' C_{60} (no ions)', ' PM6:Y6', 'PM6:Y6 (no ions)'}, 'FontSize', 25, 'Location', 'northwest')
end 

%%
if run_C60 == 1 && light == 0
    figure('Name', 'JV Plots', 'Position', [50 50 1000 1000])

    box on    
    plot(Vapp, abs(J_C60.tot(:,1)), 'LineWidth', 2, 'Color', [0 0.4470 0.7410])
    hold on
    xline(0, 'LineWidth', 2, 'Color', 'black')
    plot(Vapp, abs(J_PM6Y6.tot(:,1)), 'LineWidth', 2, 'Color', [0.4660 0.6740 0.1880])
    plot(Vapp, abs(J_PM7Y6.tot(:,1)), 'LineWidth', 2, 'Color', 'red')
    plot(Vapp, abs(J_PCE12Y6.tot(:,1)), 'LineWidth', 2, 'Color', 'black')
    hold off

    set(gca, 'FontSize', 25)
    xlabel('Voltage (V)', 'FontSize', 30)
    ylabel('Current Density (mA cm^{-2})', 'FontSize', 30)
    xlim([-0.45, 1.2])
    ylim([1e-12, 1])
    legend({'C_{60}', '', ' PM6:Y6', ' PM7:Y6', ' PCE12:Y6'}, 'FontSize', 25, 'Location', 'northwest')
end 


