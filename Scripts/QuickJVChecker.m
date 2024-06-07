% par=pc('Input_files/SnO2_MAPI_Spiro.csv');
%par=pc('Input_files/SnO2_C60_MAPI_Spiro.csv');
%par=pc('Input_files/SnO2_C60_MAPI_Spiro.csv');
%par=pc('Input_files/TiO2_MAPI_Spiro.csv');
%par=pc('Input_files/TiO2_MAPI_Spiro_TestSaP_PaperParams.csv');
% par=pc('Input_files/TiO2_MAPI_Spiro_TestSaP_4.csv');
%par=pc('Input_files/TiO2_MAPI_Spiro.csv');
% par=pc('Input_files/NiO-FACs-Al2O3-C60-Charlie.csv');
% par=pc('Input_files/EnergyOffsetSweepParameters_v5_doped.csv');
par=pc('Input_files/NiO-TripleCat-C60-Fiddled.csv');
% par = pc('Input_files/PTAA_MAPI_NegOffset_lowerVbi.csv');
% par=pc('Input_files/1_layer_test.csv');
% doped = 1;
% if doped == 1
%     par=pc('Input_files/EnergyOffsetSweepParameters_v5_doped.csv');
% elseif doped == 0
%     par=pc('Input_files/EnergyOffsetSweepParameters_v5_undoped.csv');
% end

Fiddle_with_Energetics = 0;
Fiddle_with_IonConc = 0;
IonConc = 1e18;
%%
if Fiddle_with_Energetics == 1

    %row
    DHOMO = 0.0;
    %DHOMO = Delta_HOMO(4);
    %column
    DLUMO = -0.15;
    %DLUMO = Delta_LUMO(11);
        if doped == 0
            %HTL Energetics
            par.Phi_left = -5.15;
            par.Phi_IP(1) = par.Phi_IP(3) + DHOMO;
            par.Phi_EA(1) = par.Phi_IP(1) + 2.5;
            par.EF0(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
            par.Et(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
            if par.Phi_left < par.Phi_IP(1) + 0.1
                par.Phi_left = par.Phi_IP(1) + 0.1;
            end

            %ETL Energetics
            par.Phi_right = -4.05;
            par.Phi_EA(5) = par.Phi_EA(3) + DLUMO;
            par.Phi_IP(5) = par.Phi_EA(5) - 2.5;
            par.EF0(5) = (par.Phi_IP(5)+par.Phi_EA(5))/2;
            par.Et(5) = (par.Phi_IP(5)+par.Phi_EA(5))/2;
            if par.Phi_right > par.Phi_EA(5) - 0.1
                par.Phi_right = par.Phi_EA(5) - 0.1;
            end
        
        elseif doped == 1
            %HTL Energetics
            par.Phi_left = -5.15;
            par.Phi_IP(1) = par.Phi_IP(3) + DHOMO;
            par.Phi_EA(1) = par.Phi_IP(1) + 2.5;
            par.EF0(1) = par.Phi_IP(1) + 0.1;
            par.Et(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
            if par.Phi_left < par.Phi_IP(1) + 0.1
                par.Phi_left = par.Phi_IP(1) + 0.1;
            end

            %ETL Energetics
            %Need to use opposite sign at ETL to keep energy offsets symmetric
            par.Phi_right = -4.05;
            par.Phi_EA(5) = par.Phi_EA(3) + DLUMO;
            par.Phi_IP(5) = par.Phi_EA(5) - 2.5;
            par.EF0(5) = par.Phi_EA(5) - 0.1;
            par.Et(5) = (par.Phi_IP(5) + par.Phi_EA(5))/2;
            if par.Phi_right > par.Phi_EA(5) - 0.1
                par.Phi_right = par.Phi_EA(5) - 0.1;
            end

        end
    
    par = refresh_device(par);

end

if Fiddle_with_IonConc == 1

   par.Ncat(:) = IonConc;
   par.Nani(:) = IonConc;

   par = refresh_device(par);

end

par.vsr_mode = 1;
par.frac_vsr_zone = 0.05;
par = refresh_device(par);
eqm_QJV = equilibrate(par);

%%
suns = 1;
% V_bias = 1.2;
% V_max = 1.2;
% V_min = 0;
% scan_rate = 0.2;
% deltaV = V_max - V_min;
% tmax = deltaV/scan_rate;
% 
% biased_eqm_ion = genVappStructs(eqm_QJV.ion, V_bias, 1);
% biased_eqm_el = genVappStructs(eqm_QJV.el, V_bias, 1);
% 
% illuminated_sol_ion = changeLight(biased_eqm_ion, suns, 0, 1);
% illuminated_sol_el = changeLight(biased_eqm_el, suns, 0, 1);
% 
% JV_sol_ion_rev = VappFunction(illuminated_sol_ion, 'sweep', [V_max, V_min, tmax], tmax, 200*(V_max-V_min)+1, 0);
% JV_sol_ion_fw = VappFunction(JV_sol_ion_rev, 'sweep', [V_min, V_max, tmax], tmax, 200*(V_max-V_min)+1, 0);
% 
% JV_sol_el = VappFunction(illuminated_sol_el, 'sweep', [V_max, V_min, tmax], tmax, 200*(V_max-V_min)+1, 0);

JV_sol_ion = doCV(eqm_QJV.ion, suns, -0.2, 1.3, -0.2, 1e-4, 1, 301);
%JV_sol_el = doCV(eqm_QJV.el, suns, -0.2, 1.3, -0.2, 0.1, 1, 301);
% JV_sol_ion = doCV(eqm_QJV.ion, suns, -0.2, 1.2, -0.2, 100e-3, 1, 281);

Plot_Current_Contributions(JV_sol_ion) 
%Plot_Current_Contributions(JV_sol_el) 
stats_ion = CVstats(JV_sol_ion)
%stats_el = CVstats(JV_sol_el)

% %% Plot JVs
% figure('Name', 'JVPlot', 'Position', [100 100 1250 1250])
% % colors_JV = {[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],[0 0.4470 0.7410],[0.9290 0.6940 0.1250]};
% % % v_fw = dfana.calcVapp(JV_sol_ion_fw);
% % % v_rev = dfana.calcVapp(JV_sol_ion_rev);
% v = dfana.calcVapp(JV_sol_el)';
% % v_ion = dfana.calcVapp(JV_sol_ion);
% % 
% % hold on
% % xline(0, 'black', 'LineWidth', 1)
% % yline(0, 'black', 'LineWidth', 1)
% % 
% % % j_fw = dfana.calcJ(JV_sol_ion_fw).tot(:,1);
% % % j_rev = dfana.calcJ(JV_sol_ion_rev).tot(:,1);
% j_el = dfana.calcJ(JV_sol_el).tot(:,1);
% j_ion = dfana.calcJ(JV_sol_ion).tot(:,1);
% % plot(v_fw(1:end), j_fw(1:end)*1000, 'color', [0.4660 0.6740 0.1880], 'LineWidth', 3) 
% plot(v_ion(1:end), j_ion(1:end)*1000, 'color', [0.4660 0.6740 0.1880], 'LineWidth', 3) 
% hold on
% % plot(v_rev(1:end), j_rev(1:end)*1000, 'color', [0.4660 0.6740 0.1880], 'LineWidth', 3) 
% plot(v_el(1:end), j_el(1:end)*1000, 'color', [0 0.4470 0.7410], 'LineWidth', 3)
% hold on
% 
% hold off
% 
% box on 
% set(gca, 'FontSize', 25)
% xlim([-0.15, 1.2])
% ylim([-25,5])
% legend({'','','Mobile Ions','No Mobile Ions'}, 'Location', 'northwest', 'FontSize', 30)
% xlabel('Voltage(V)', 'FontSize', 30)
% ylabel('Current Density (mAcm^{-2})', 'FontSize', 30)
% ax1 = gcf;

%% Save Plots at 300 dpi
save_plot = 0;

if save_plot == 1 
    exportgraphics(ax1, ...
    'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\LSR\0p15ev_Symmetric_Offest_JV.png', ...
    'Resolution', 300)
end 

%%
%Make one sun solution at a given applied voltage
run = 0; 
if run == 1
    %Vapp = 0.62; %Uniform ion distribution from JV for negtive offset case
    Vapp = stats_ion.Voc_f; %Uniform ion distribution from JV
    sol_ill = changeLight(eqm_QJV.ion, 1.1, 0, 1);
    sol_ill_bias = genVappStructs(sol_ill, Vapp, 1);
    try
        CV_sol_startbias = doCV(sol_ill_bias, 1.1, Vapp, -0.20, 1.20, 100, 1, 281);
    catch
        warning('No joy.')
    end

    Plot_Current_Contributions(CV_sol_startbias)
    stats_ions_bias = CVstats(CV_sol_startbias)

    dfplot.ELxnpxacx(JV_sol_ion, 1e4*(0.2+Vapp))
    dfplot.ELxnpxacx(sol_ill_bias, sol_ill_bias.t(end))
end     


