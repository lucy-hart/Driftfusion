%% File to test that doSaP is doing what I intend
par = pc('Input_files/EnergyOffsetSweepParameters_v5_doped_asymmetric.csv');
% par=pc('Input_files/SnO2_MAPI_Spiro_TestSaP.csv');
% par=pc('Input_files/TiO2_MAPI_Spiro_TestSaP.csv');
% par = pc('Input_files/PTAA_MAPI_NegOffset_lowerVbi.csv');

compare_fixed_ion_JV = 1;

DHOMO = 0.0001;
DLUMO = -0.15;
           
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
par.Phi_right = -4.05;
par.Phi_EA(5) = par.Phi_EA(3) + DLUMO;
par.Phi_IP(5) = par.Phi_EA(5) - 2.5;
par.EF0(5) = par.Phi_EA(5) - 0.1;
par.Et(5) = (par.Phi_IP(5) + par.Phi_EA(5))/2;
if par.Phi_right > par.Phi_EA(5) - 0.1
    par.Phi_right = par.Phi_EA(5) - 0.1;
end

% par.frac_vsr_zone = 0.05;
par.RelTol_vsr = 0.1;
par = refresh_device(par);

eqm = equilibrate(par);

%% See what device performance is at illumination used for SaP measurement
check_JV = 1;
suns = 1;
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
%Vbias = linspace(0,1.2,13);
Vbias = [1.18];
Vpulse = linspace(0,1.2,49);
tramp = 8e-4;
tsample = 1e-3;
tstab = 120;

sol = doSaP_v2(eqm.ion, Vbias, Vpulse, tramp, tsample, tstab, suns);

%% Do JVs with mobseti = 0 to compare the SaP JVs
fixed_ion_JVs = cell(length(Vbias));
J_fixed_ion = cell(length(Vbias));
if compare_fixed_ion_JV == 1
    for i=1:length(Vbias)
        disp(['Doing JV for Vstab = ' num2str(Vbias(i)) ' V'])
        sol{i,1}.par.mobseti = 0;
        fixed_ion_JVs{i} = doCV(sol{i,1}, suns, -0.2, 1.2, -0.2, 1e-3, 1, 281);
        J_fixed_ion{i} = dfana.calcJ(fixed_ion_JVs{i}).tot(:,1);
        if i == 1
            V_fixed_ion = dfana.calcVapp(fixed_ion_JVs{i});
        end
    end
end

%% Look at the voltage stabilisation
t = sol{1, 1}.t;

Jstab = zeros(length(Vbias), length(t));

for i = 1:length(Vbias)
    Jtemp = dfana.calcJ(sol{i,1});
    Jstab(i,:) = Jtemp.tot(:,1);
end

%%
figure('Name', 'JstabData')
hold on
box on
for i = 1:length(Vbias)
    plot(t, 1e3*Jstab(i,:), 'DisplayName', [num2str(Vbias(i)) ' V'])
end
ylabel('Current Density (mA cm^{-2})')
xlabel('Time (s)')
legend()
%% Extract the current values from the pulsed JVs
bias = 1;
t = sol{bias, 2}.t;

Jt = zeros(length(Vpulse), length(t));

for i = 1:length(Vpulse)
    Jtemp = dfana.calcJ(sol{bias,i+1});
    Jtot = Jtemp.tot(:,1);
    if length(Jtot) < length(t)
        Jtot(end+1:numel(t)) = Jtot(1);
    end
    Jt(i,:) = Jtot(:,1);
end

%%
figure('Name', 'PulsedJVCurrents')
hold on
box on
for i = 1:length(Vpulse)
    plot(t, 1e3*Jt(i,:), 'DisplayName', [num2str(Vpulse(i)) ' V'])
end
ylabel('Current Density (mA cm^{-2})')
xlabel('Time (s)')
legend()

%% Plot pulsed JVs 
figure('Name', 'PulsedJVs')
cmap = colormap(parula(length(Vbias)));
hold on
box on
xline(0, 'black', 'HandleVisibility', 'off')
yline(0, 'black', 'HandleVisibility', 'off')
for i = 1:length(Vbias)
    Jpulse = zeros(1, length(Vpulse));
    Jpulse2 = zeros(1, length(Vpulse));
    for j = 1:length(Vpulse)
        Jpulse(j) = sol{i,j+1}.Jpulse;
    end
    plot(Vpulse(:), 1e3*Jpulse, 'DisplayName', num2str(Vbias(i), '%.1f'), 'color', cmap(i,:))
    if compare_fixed_ion_JV == 1
        plot(V_fixed_ion, 1e3*J_fixed_ion{i}, 'HandleVisibility', 'Off', 'color', 'black', 'LineStyle', ':')
    end
end
plot(v(1:141), J_el(1:141)*1000, 'color', 'black', 'LineWidth', 3, 'LineStyle', ':')
ylabel('Current Density (mA cm^{-2})')
ylim([-25, 10])
xlabel('Voltage (V)')
xlim([Vpulse(1), 1.2])
legend()
title(legend, 'V_{bias} (V)')

%%
i = 1;
Jpulse = zeros(length(Vpulse),1);
for j = 1:length(Vpulse)
    Jpulse(j) = sol{i,j+1}.Jpulse;
end

