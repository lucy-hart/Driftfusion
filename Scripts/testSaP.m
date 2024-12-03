%% File to test that doSaP is doing what I intend
% par = pc('Input_files/EnergyOffsetSweepParameters_v5_doped.csv');
%par = pc('Input_files/EnergyOffsetSweepParameters_v5_undoped.csv');
%par=pc('Input_files/SnO2_MAPI_Spiro_TestSaP.csv');
par=pc('Input_files/TiO2_MAPI_Spiro_TestSaP_4.csv');
% par=pc('Input_files/NiO-TripleCat-C60.csv');
% par=pc('Input_files/TiO2_MAPI_Spiro_TestSaP_PaperParams.csv');
% par.prob_distro_function = 'Boltz';
% par = pc('Input_files/PTAA_MAPI_NoOffset.csv');
par.RelTol_vsr = 0.1;
compare_fixed_ion_JV = 1;

% DHOMO = 0.35;
% DLUMO = -0.35;
% IonConc = 1e17;
%            
% %HTL Energetics
% par.Phi_left = -5.15; 
% par.Phi_IP(1) = par.Phi_IP(3) + DHOMO;
% par.Phi_EA(1) = par.Phi_IP(1) + 2.5;
% par.EF0(1) = par.Phi_IP(1) + 0.1;
% par.Et(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
% if par.Phi_left < par.Phi_IP(1) + 0.1
%     par.Phi_left = par.Phi_IP(1) + 0.1;
% end
% %ETL Energetics
% par.Phi_right = -4.05;
% par.Phi_EA(5) = par.Phi_EA(3) + DLUMO;
% par.Phi_IP(5) = par.Phi_EA(5) - 2.5;
% par.EF0(5) = par.Phi_EA(5) - 0.1;
% par.Et(5) = (par.Phi_IP(5) + par.Phi_EA(5))/2;
% if par.Phi_right > par.Phi_EA(5) - 0.1
%     par.Phi_right = par.Phi_EA(5) - 0.1;
% end
% 
% par.Ncat(:) = IonConc;
% par.Nani(:) = IonConc;
% 
% % par.frac_vsr_zone = 0.05;
% par.RelTol_vsr = 0.1;
% par = refresh_device(par);

eqm = equilibrate(par);

%% See what device performance is at illumination used for SaP measurement
check_JV = 0;
suns = 0;
if check_JV ==1 
    JVsol_el = doCV(eqm.el, suns, -0.2, 1.2, -0.2, 1e-4, 1, 281);
    JVsol_ion = doCV(eqm.ion, suns, -0.2, 1.2, -0.2, 1e-4, 1, 281);
    
    figure('Name', 'JVPlot')
    v = dfana.calcVapp(JVsol_el);
    J_el = dfana.calcJ(JVsol_el).tot(:,1);
    J_ion = dfana.calcJ(JVsol_ion).tot(:,1);
    box on
    hold on
    xline(0, 'color', 'black', 'LineWidth', 1)
    yline(0, 'color', 'black', 'LineWidth', 1)
    plot(v(:), J_el(:)*1000, 'color', [0.4660 0.6740 0.1880], 'LineWidth', 3)
    plot(v(:), J_ion(:)*1000, 'color', 'black', 'LineWidth', 3)
    hold off
    xlim([-0.2, 1.2])
    xlabel('Current Density (mA cm^{-2})')
    ylim([-25, 5])
    ylabel('Voltage (V)')
end
%% Do the SaP measurement
Vbias = linspace(0,1.5,16);
%Vbias = [0 0.1];
%Vpulse = linspace(0,1.3,27);
Vpulse = [0];
tramp = 8e-4;
tsample = 1e-3;
tstab = 200;

sol = doSaP_v2(eqm.ion, Vbias, Vpulse, tramp, tsample, tstab, suns, 1);

%% Do JVs with mobseti = 0 to compare the aP JVs
fixed_ion_JVs = cell(1, length(Vbias));
J_fixed_ion = cell(1, length(Vbias));
if compare_fixed_ion_JV == 1
    for i=1:length(Vbias)
        disp(['Doing JV for Vstab = ' num2str(Vbias(i)) ' V'])
        sol{i,1}.par.mobseti = 0;
        try
            fixed_ion_JVs{i} = doCV(sol{i,1}, suns, -0.1, 1.05, -0.1, 1, 0.5, 116);
            J_fixed_ion{i} = dfana.calcJ(fixed_ion_JVs{i}).tot(:,1);
        catch
            warning(['Fixed ion JV failed at Vbias of ', num2str(Vbias(i))])
            J_fixed_ion{i} = zeros(116,1);
        end
        if i == 1
            V_fixed_ion = dfana.calcVapp(fixed_ion_JVs{i});
        end
    end
end

%% Look at the voltage stabilisation
% t = sol{1, 1}.t;
% 
% Jstab = zeros(length(Vbias), length(t));
% 
% for i = 1:length(Vbias)
%     Jtemp = dfana.calcJ(sol{i,1});
%     Jstab(i,:) = Jtemp.tot(:,1);
% end

%%
% figure('Name', 'JstabData')
% hold on
% box on
% for i = 1:length(Vbias)
%     plot(t, 1e3*Jstab(i,:), 'DisplayName', [num2str(Vbias(i)) ' V'])
% end
% ylabel('Current Density (mA cm^{-2})')
% xlabel('Time (s)')
% legend()
%% Extract the current values from the pulsed JVs
bias = 1;
t = sol{bias, 2}.t;

Jt = zeros(length(Vpulse), length(t));

for i = 1:length(Vpulse)
    try
        Jtemp = dfana.calcJ(sol{bias,i+1});
        Jtot = Jtemp.tot(:,1);
    catch
        Jtot = 0;
    end
    if length(Jtot) < length(t)
        Jtot(end+1:numel(t)) = Jtot(1);
    end
    Jt(i,:) = Jtot(:,1);
end

%%
% figure('Name', 'PulsedJVCurrents')
% hold on
% box on
% for i = 1:length(Vpulse)
%     plot(t, 1e3*Jt(i,:), 'DisplayName', [num2str(Vpulse(i)) ' V'])
% end
% ylabel('Current Density (mA cm^{-2})')
% xlabel('Time (s)')
% legend()

%% Plot pulsed JVs 
figure('Name', 'PulsedJVs')
cmap = colormap(parula(length(Vbias)));
cmap = flip(cmap);
hold on
box on
xline(0, 'black', 'HandleVisibility', 'off')
yline(0, 'black', 'HandleVisibility', 'off')

for i = 1:length(Vbias)
    Jpulse = zeros(1, length(Vpulse));
    for j = 1:length(Vpulse)
        Jpulse(j) = sol{i,j+1}.Jpulse;
    end
    plot(Vpulse(Jpulse ~= 0), 1e3*Jpulse(Jpulse ~= 0), 'DisplayName', num2str(Vbias(i), '%.2f'), 'color', cmap(i,:))
    if compare_fixed_ion_JV == 1
        plot(V_fixed_ion, 1e3*J_fixed_ion{i}, 'HandleVisibility', 'Off', 'color', cmap(i,:), 'LineStyle', '-')
    end
end

if check_JV == 1
    plot(v(1:141), J_el(1:141)*1000, 'color', 'black', 'LineWidth', 3, 'LineStyle', '-')
end
ylabel('Current Density (mA cm^{-2})')
ylim([-25, 5])
xlabel('Voltage (V)')
xlim([0, 1.05])
legend()
%title(legend, 'V_{bias} (V)')

%% Plot Vx at 0.6 V
% figure('Name', 'Vx-Vbias')
% cmap = colormap(parula(length(Vbias)));
% hold on
% box on
% xline(0, 'black', 'HandleVisibility', 'off')
% yline(0, 'black', 'HandleVisibility', 'off')
% 
% for i = 1:length(Vbias)
%     plot((fixed_ion_JVs{i}.x*1e7), flip(fixed_ion_JVs{i}.u(81,:,1)), ...
%         'DisplayName', num2str(Vbias(i), '%.1f'), 'color', cmap(i,:))
% end
% 
% xlabel('Position (nm)')
% ylim([-0.1, 0.5])
% ylabel('Electrostaic Potential (V)')
% xlim([0, 1e7*max(fixed_ion_JVs{1}.x)])
% legend()
% title(legend, 'V_{bias} (V)')
% 
%% Do SaP analysis
% all_data = {J_fixed_ion};%, sol2};
% 
% idx_stop = floor(length(V_fixed_ion));
% V_SaP_Analysis = V_fixed_ion(1:idx_stop);
% Jvalues = zeros(length(Vbias), idx_stop, length(all_data));
% Voc = zeros(length(all_data), length(Vbias));
% dJdV_Voc = zeros(length(all_data), length(Vbias));
% for k = 1:length(all_data)
%     data = all_data{k};
%     for i = 1:length(Vbias)
%         %for j = 1:length(V_fixed_ion)
%             Jvalues(i,:,k) = data{i}(1:idx_stop);
%         %end
%     end
% end  
% % 
% colours = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880], [0.9290 0.6940 0.1250]};
% h = Vbias(2) - Vbias(1);
% for k = 1:length(all_data)
%     for i = 1:length(Vbias)
%         J_temp = Jvalues(i,:,k);
%         Voc(i) = interp1(J_temp(J_temp~=0), V_SaP_Analysis(J_temp~=0), 0);
%         dJdV = gradient(J_temp(J_temp~=0), V_SaP_Analysis(J_temp~=0));
% %         plot(V_SaP_Analysis(J_temp~=0), dJdV, color = colours{k})
%         hold on
%         dJdV_Voc(k,i) = interp1(V_SaP_Analysis(J_temp~=0), dJdV, Voc(i));
%     end
% end

%% Plot SaP analysis
%NB: Smaller DeltaE_ETL also had vs = 50 cms-1
% figure('Name', 'SaP-Analysis')
% 
% subplot(1,1,1)
% hold on
% box on
% xline(0, 'black', 'HandleVisibility', 'off')
% yline(0, 'black', 'HandleVisibility', 'off')
% %plot(Vbias, dJdV_Voc(2,:),'DisplayName','v_s = 0.05 cm s^{-1}', 'Color', [0 0.4470 0.7410])
% plot(Vbias, dJdV_Voc(1,:),'DisplayName','v_s = 50 cm s^{-1}', 'Color', [0.8500 0.3250 0.0980])
% 
% xlabel('V_{bias} (V)')
% xlim([Vbias(1), Vbias(end)])
% ylabel('dJ/dV|_{Voc} (\Omega^{-1} cm^{-2})')
% legend('Location', 'northwest')

% subplot(1,2,2)
% hold on
% box on
% xline(0, 'black', 'HandleVisibility', 'off')
% yline(0, 'black', 'HandleVisibility', 'off')
% plot(Vbias, dJdV_Voc(2,:),'DisplayName','V_{flat} = 0.6 eV', 'Color', [0 0.4470 0.7410])
% plot(Vbias, dJdV_Voc(3,:),'DisplayName','V_{flat} = 0.8 eV (asymmetric)', 'Color', [0.4660 0.6740 0.1880])
% plot(Vbias, dJdV_Voc(4,:),'DisplayName','V_{flat} = 0.8 eV (symmetric)', 'Color', [0.9290 0.6940 0.1250])
% 
% xlabel('V_{bias} (V)')
% xlim([Vbias(1), Vbias(end)])
% ylabel('dJ/dV|_{Voc} (\Omega^{-1} cm^{-2})')
% legend('Location', 'northwest')

