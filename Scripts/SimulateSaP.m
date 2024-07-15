%% File to Simulate SaP Measurements
TripleCat = 1;

if TripleCat == 1
    par=pc('Input_files/NiO-TripleCat-C60.csv');
    suns = 0.65;
    Vmax = 1.3;
    num_points = 281;
else
    par=pc('Input_files/TiO2_MAPI_Spiro.csv');
    suns = 1;
    Vmax = 1.1;
    num_points = 241;
end 
par.RelTol_vsr = 0.1;
compare_fixed_ion_JV = 1;

eqm = equilibrate(par);

%% Do the SaP measurement
Vbias = linspace(0,1.3,14);
Vpulse = linspace(0,1.3,27);
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
            fixed_ion_JVs{i} = doCV(sol{i,1}, suns, -0.1, Vmax, -0.1, 1, 0.5, num_points);
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
        plot(V_fixed_ion, 1e3*J_fixed_ion{i}, 'HandleVisibility', 'Off', 'color', 'Black', 'LineStyle', '--')
    end
end

ylabel('Current Density (mA cm^{-2})')
ylim([-25, 5])
xlabel('Voltage (V)')
xlim([0, 1.05])
legend()
title(legend, 'V_{bias} (V)')
