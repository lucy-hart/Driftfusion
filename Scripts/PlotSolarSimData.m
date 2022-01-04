%Plot JV curves measured from the SolarSim

%% Read in data (and get rid of noisy points)
path = 'C:\Users\User\OneDrive - Imperial College London\PhD\Figures\Weidong_ETL\SolarSimData';
Kloc6 = readtable(join([path,'\Kloc6\fw_102_6'], ''));
outlier = Kloc6.Voltage == 0.2;
Kloc6(outlier, :) = [];
PCBM = readtable(join([path,'\PCBM\788_fw_100_1'], ''));
outlier = PCBM.Voltage == 0.86;
PCBM(outlier, :) = [];
ICBA = readtable(join([path,'\ICBA\fw_100_1'], ''));

%% Do plotting
figure(42)
plot(Kloc6{:,1}, Kloc6{:,3}, 'r', PCBM{:,1}, PCBM{:,3}, 'g', ICBA{:,1}, ICBA{:,3}, 'b')
legend('Kloc-6', 'PCBM', 'ICBA', 'Location', 'northwest')
xlabel('Voltage (V)')
ylabel('Current Density (mAcm^{-2}')
xlim([-0, 1.3])
ylim([-30, 0])
