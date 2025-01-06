%% File to test how ion mobility affects accuracy of SaP method 
par=pc('Input_files/NiO-TripleCat-C60.csv');
par.RelTol_vsr = 0.08;

%% Do the SaP measurement
Vbias = linspace(0,1.2,13);
Vpulse = linspace(0,1.3,14);
%Vpulse = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.0 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0];
%mobilities = logspace(-12,-9,4);
mobilities = [1e-10 1e-9 1e-8];
%mobilities = [1e-9];
tramp = 8e-4;
tsample = 1e-3-tramp;
tstab = 120;
sol_pulsed =cell(1,length(mobilities));
suns=1;

for i = 1:length(mobilities)
    idx_active = par.active_layer;
    par.mu_c(idx_active) = mobilities(i);
    par.mu_a(idx_active) = mobilities(i);
    par = refresh_device(par);
    eqm = equilibrate(par);
    sol_pulsed{i} = doSaP_v3(eqm.ion, Vbias, Vpulse, tramp, tsample, tstab, suns, 1, 1);
end

%% Do JVs with mobseti = 0 to compare the SaP JVs
fixed_ion_JVs = cell(length(mobilities), length(Vbias));
J_fixed_ion = cell(length(mobilities), length(Vbias));
compare_fixed_ion_JV = 1;

if compare_fixed_ion_JV == 1
    for j = 1:length(mobilities)
        sol = sol_pulsed{j};
        for i=1:length(Vbias)
            disp(['Doing JV for Vstab = ' num2str(Vbias(i)) ' V'])        
            sol{i,1}.par.mobseti = 0;
            try
                fixed_ion_JVs{j,i} = doCV(sol{i,1}, suns, -0.1, max(Vpulse), -0.1, 1, 0.5, 1+int32(2*(100*(max(Vpulse)+0.1))));
                J_fixed_ion{j,i} = dfana.calcJ(fixed_ion_JVs{j,i}).tot(:,1);
            catch
                warning(['Fixed ion JV failed at Vbias of ', num2str(Vbias(i))])
                J_fixed_ion{j,i} = zeros(1+int32(2*(10*(Vpulse(end)+0.1))),1);
            end
            if i == 1 && j == 1
                V_fixed_ion = dfana.calcVapp(fixed_ion_JVs{i});
            end
        end
    end
end

%% Look at the voltage stabilisation
which = 3;
sol = sol_pulsed{which};
t = sol{1, 1}.t;

Jstab = zeros(length(Vbias), length(t));

for i = 1:length(Vbias)
    Jtemp = dfana.calcJ(sol{i,1});
    Jstab(i,:) = Jtemp.tot(:,1);
end
figure('Name', 'JstabData')
cmap = colormap(parula(length(Vbias)));
cmap = flip(cmap);

subplot(1,2,1)
hold on
box on
for i = 1:length(Vbias)
    plot(t, 1e3*Jstab(i,:), 'DisplayName', [num2str(Vbias(i)) ' V'], 'color', cmap(i,:))
end
ylabel('Current Density (mA cm^{-2})')
xlabel('Time (s)')
legend()

subplot(1,2,2)
hold on
box on
xline(0, 'black', 'HandleVisibility', 'off')
yline(0, 'black', 'HandleVisibility', 'off')

for i = 1:length(Vbias)
    Jpulse = zeros(1, length(Vpulse));
    for j = 1:length(Vpulse)
        Jpulse(j) = sol{i,j+1}.Jpulse;
    end
    plot(Vpulse(Jpulse ~= 0), 1e3*Jpulse(Jpulse ~= 0), 'DisplayName', num2str(Vbias(i), '%.2f'), 'color', cmap(i,:), 'LineStyle', '-', 'Marker','o')
%     if compare_fixed_ion_JV == 1
%         plot(V_fixed_ion, 1e3*J_fixed_ion{which,i}, 'HandleVisibility', 'Off', 'color', 'red', 'LineStyle', '-')
%     end
end

ylabel('Current Density (mA cm^{-2})')
ylim([-25, 5])
xlabel('Voltage (V)')
xlim([0, 1.3])
legend()
%title(legend, 'V_{bias} (V)')



