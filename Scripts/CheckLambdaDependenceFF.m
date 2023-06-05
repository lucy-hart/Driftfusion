parC60 = pc('Input_files/HTL_MAPI_C60.csv');
parPM6Y6 = pc('Input_files/HTL_MAPI_PM6Y6.csv');

devices = {parC60, parPM6Y6};
num_devices = length(devices);
wavelengths = [404, 750];
power = [54, 17];
eqm = cell(2,2);
CV_sols = cell(2,2);
CV_stats = cell(2,2);
JValues = cell(2,2);

for j = 1:length(wavelengths)
    for i = 1:num_devices    
        devices{i}.light_source1 = 'laser';
        devices{i}.laser_lambda1 = wavelengths(j);
        devices{i}.pulsepow = power(j);
        devices{i} = refresh_device(devices{i});
        eqm{i,j} = equilibrate(devices{i});
        CV_sols{i,j} = doCV(eqm{i,j}.ion, 1.0, -0.2, 1.2, -0.2, 10e-3, 1, 281);
        CV_stats{i,j} = CVstats(CV_sols{i,j});
        JValues{i,j} = dfana.calcJ(CV_sols{i,j});
    end
end

Vapp = dfana.calcVapp(CV_sols{1,1});

%%
figure('Name', 'JV Plots', 'Position', [50 50 1000 1000])
norm = 0;

for k = 1:num_devices
    subplot(1,2,k)
    box on
    hold on
    xline(0, 'LineWidth', 2, 'Color', 'black')
    yline(0, 'LineWidth', 2, 'Color', 'black')
    if norm == 0
        plot(Vapp, 1e3*JValues{k,1}.tot(:,1), 'LineWidth', 4, 'Color', 'b')
        plot(Vapp, 1e3*JValues{k,2}.tot(:,1), 'LineWidth', 4, 'Color', 'r')
    elseif norm == 1
        plot(Vapp, -JValues{k,1}.tot(:,1)/min(JValues{k,1}.tot(:,1)), 'LineWidth', 4, 'Color', 'b')
        plot(Vapp, -JValues{k,2}.tot(:,1)./min(JValues{k,2}.tot(:,1)), 'LineWidth', 4, 'Color', 'r')    
    end
    set(gca, 'FontSize', 25)
    xlabel('Voltage (V)', 'FontSize', 30)
    xlim([-0.15, 1.2])
    if norm == 0
        ylabel('Current Density (mA cm^{-2})', 'FontSize', 30)
        ylim([-20, 5])
    elseif norm == 1
        ylabel('J/J_{SC}', 'FontSize', 30)
        ylim([-1.1, 0.2])
    end
    legend({'', '', '404 nm', '750 nm'}, 'FontSize', 25, 'Location', 'northwest')
end
    hold off


