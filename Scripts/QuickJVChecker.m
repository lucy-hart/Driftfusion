%par=pc('Input_files/EnergyOffsetSweepParameters_v2.csv');
%par=pc('Input_files/EnergyOffsetSweepParameters_v3.csv');
par=pc('Input_files/SnO2_MAPI_Spiro.csv');
%par = pc('Input_files/PTAA_MAPI_NegOffset.csv');
%par = pc('Input_files/PTAA_MAPI_NegOffset_lowerVbi.csv');
%par = pc('Input_files/PTAA_MAPI_NoOffset.csv');
%par = pc('Input_files/PTAA_MAPI_PosOffset.csv');
%par = pc('Input_files/PTAA_MAPI_PCBM_ForPaper.csv');

doped = 5;
if doped == 1
    par=pc('Input_files/EnergyOffsetSweepParameters_v4_doped.csv');
elseif doped == 0
    par=pc('Input_files/EnergyOffsetSweepParameters_v4_undoped.csv');
end

Fiddle_with_Energetics = 0;
Fiddle_with_IonConc = 0;
IonConc = 1e19;
%%
if Fiddle_with_Energetics == 1

    %row
    DHOMO = 0.3;
    %DHOMO = Delta_HOMO(4);
    %column
    DLUMO = -0.3;
    %DLUMO = Delta_LUMO(11);
        if doped == 0
            %HTL Energetics
            par.Phi_left = -5.15;
            par.Phi_IP(1) = par.Phi_IP(3) + DHOMO;
            par.Phi_EA(1) = par.Phi_IP(1) + 2.5;
            par.EF0(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
            par.Et(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
            if par.Phi_left < par.Phi_IP(1) + 0.01
                par.Phi_left = par.Phi_IP(1) + 0.01;
            end

            %ETL Energetics
            par.Phi_right = -4.05;
            par.Phi_EA(5) = par.Phi_EA(3) + DLUMO;
            par.Phi_IP(5) = par.Phi_EA(5) - 2.5;
            par.EF0(5) = (par.Phi_IP(5)+par.Phi_EA(5))/2;
            par.Et(5) = (par.Phi_IP(5)+par.Phi_EA(5))/2;
            if par.Phi_right > par.Phi_EA(5) - 0.01
                par.Phi_right = par.Phi_EA(5) - 0.01;
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

eqm_QJV = equilibrate(par);

%%
suns = 1;

illuminated_sol_ion = changeLight(eqm_QJV.ion, suns, 0, 1);
illuminated_sol_el = changeLight(eqm_QJV.el, suns, 0, 1);

V_bias = 1.5;
V_max = 1.1;
V_min = 0;

biased_eqm_ion = genVappStructs(illuminated_sol_ion, V_bias, 1);
biased_eqm_el = genVappStructs(illuminated_sol_el, V_bias, 1);

JV_sol_ion = doCV(biased_eqm_ion, suns, V_bias, V_max, V_min, 1e-3, 1, 221);
JV_sol_el = doCV(biased_eqm_el, suns, V_bias, V_max, V_min, 1e-3, 1, 221);

Plot_Current_Contributions(JV_sol_ion) 
Plot_Current_Contributions(JV_sol_el) 
stats_ion = CVstats(JV_sol_ion)
stats_el = CVstats(JV_sol_el)

%% Plot JVs
figure('Name', 'JVPlot', 'Position', [100 100 1250 1250])
colors_JV = {[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],[0 0.4470 0.7410],[0.9290 0.6940 0.1250]};
v = dfana.calcVapp(JV_sol_ion);
v_el = dfana.calcVapp(JV_sol_el);

hold on
plot(v(:), zeros(1,length(v)), 'black', 'LineWidth', 1)
plot(zeros(1,50), linspace(-30, 10, 50), 'black', 'LineWidth', 1)

j = dfana.calcJ(JV_sol_ion).tot(:,1);
j_el = dfana.calcJ(JV_sol_el).tot(:,1);
plot(v(:), j(:)*1000, 'color', [0.4660 0.6740 0.1880], 'LineWidth', 3) 
hold on
plot(v_el(1:116), j_el(1:116)*1000, '--', 'color', [0.4660 0.6740 0.1880], 'LineWidth', 3)
hold on

hold off

box on 
set(gca, 'FontSize', 25)
xlim([-0.15, 1.2])
ylim([-25,5])
legend({'','','Mobile Ions','No Mobile Ions'}, 'Location', 'northwest', 'FontSize', 30)
xlabel('Voltage(V)', 'FontSize', 30)
ylabel('Current Density (mAcm^{-2})', 'FontSize', 30)
ax1 = gcf;

%% Save Plots at 300 dpi
save = 0;

if save == 1 
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


