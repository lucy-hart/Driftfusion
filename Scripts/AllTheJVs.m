%% Plot JVs
%Run after Weidong_ETL.m
figure('Name', 'JVPlot', 'Position', [100 100 1250 2000])
%Colour order is green, blue, red, yellow, purple
colors_allJV = {[0.4660 0.6740 0.1880],[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250],[0.4940 0.1840 0.5560]};

box on 
hold on
for m = 1:num_devices
    v = dfana.calcVapp(CV_solutions_ion{m});
    v_el = dfana.calcVapp(CV_solutions_el{m});
    j = -dfana.calcJ(CV_solutions_ion{m}).tot(:,1);
    j_el = -dfana.calcJ(CV_solutions_el{m}).tot(:,1);

    plot(v(:), j(:)*1000, 'color', colors_allJV{m}, 'LineWidth', 3) 
    plot(v_el(1:151), j_el(1:151)*1000, '--', 'color', colors_allJV{m}, 'LineWidth', 3)
end
hold off

set(gca, 'FontSize', 25)
xlim([0, 1.3])
ylim([0,27])
xlabel('Voltage (V)', 'FontSize', 30)
ylabel('Current Density (mAcm^{-2})', 'FontSize', 30)
legend({'  PS1', '', '  PS2', '',  '  PS3', '',  '  PS4', '',  '  PS5' ''}, 'Location', 'southwest', 'FontSize', 30)
ax1 = gca;

%% Save Plot
save = 0;

if save == 1 
    exportgraphics(ax1, ...
    'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Weidong_ETL\Paper\v2\AlltheJVs.png', ...
    'Resolution', 300)
end 