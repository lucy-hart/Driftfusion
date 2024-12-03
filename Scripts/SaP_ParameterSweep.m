%% File to test that what varying device parameters does to SaP measurements
par=pc('Input_files/TripleCat-ForSaP.csv');
par.RelTol_vsr = 0.1;

%IonConc = [1e16 5e16 1e17 5e17 1e18];
% vsurf = [5 50 500];
%vsurf = [1e7 1e5 1e3 10];
%tau = [10e-9 50e-9 100e-9 500e-9 1e-6];
Vbi_offset = [0 0.05 0.1 0.15];
num_samples = length(Vbi_offset);
devices = cell(num_samples, 1);

for i = 1:num_samples
    if i == 1
        Phi0_left = par.Phi_left;
        Phi0_right = par.Phi_right;
    end
%     par.Ncat(:) = IonConc(i);
%     par.Nani(:) = IonConc(i);
    par.Ncat(:) = 1e17;
    par.Nani(:) = 1e17;
%     par.sn_r = vsurf(i);
%     par.sp_r = vsurf(i);
%     par.sn_l = 1;
%     par.sp_r = 1;
%     par.sn_r = 1e7;
%     par.sp_l = 1e7;
    par.Phi_left = Phi0_left + Vbi_offset(i);
    par.Phi_right = Phi0_right - Vbi_offset(i);
%     par.taun = tau(i);
%     par.taup = tau(i);
    par = refresh_device(par);
    devices{i} = equilibrate(par);
end

%% Get the 1 Sun QSS solutions at different biases
sol = cell(num_samples, 1);
Vbias = linspace(0,1.5,16);
suns = 1;
%Vbias(15) = [];
%Vbias = linspace(0,0.2,3);
Vpulse = [0];
tramp = 8e-4;
tsample = 1e-3;
tstab = 200;

for i = 1:num_samples
    sol{i} = doSaP_v2(devices{i}.ion, Vbias, Vpulse, tramp, tsample, tstab, suns, 0);
end
%% Do JVs with mobseti = 0 
fixed_ion_JVs = cell(num_samples, length(Vbias));
J_fixed_ion = cell(num_samples, length(Vbias));
Vmax = 1.2;
Vmin = 0;
numpoints = 100*(Vmax-Vmin) + 1;

for j=1:num_samples
    for i=1:length(Vbias)
        sol_temp = sol{j}{i};
        disp(['Doing JV for Vstab = ' num2str(Vbias(i)) ' V and j = ', num2str(j)])
        sol_temp.par.mobseti = 0;
        try
            fixed_ion_JVs{j,i} = doCV(sol_temp, suns, Vmin, Vmax, Vmin, 1, 0.5, numpoints);
            J_fixed_ion{j,i} = dfana.calcJ(fixed_ion_JVs{j,i}).tot(:,1);
        catch
            warning(['Fixed ion JV failed at Vbias = ', num2str(Vbias(i)), ' V for j = ', num2str(j)])
            J_fixed_ion{j,i} = zeros(numpoints,1);
        end
        if i == 1 && j == 1
            V_fixed_ion = dfana.calcVapp(fixed_ion_JVs{j,i});
        end
    end
end

%% Plot pulsed JVs 
num = 3;
figure('Name', 'PulsedJVs')
cmap = colormap(parula(length(Vbias)));
cmap = flip(cmap);
hold on
box on
xline(0, 'black', 'HandleVisibility', 'off')
yline(0, 'black', 'HandleVisibility', 'off')

for i = 1:length(Vbias)
    plot(V_fixed_ion, 1e3*J_fixed_ion{num,i}, 'HandleVisibility', 'Off', 'color', cmap(i,:), 'LineStyle', '-')
end

ylabel('Current Density (mA cm^{-2})')
ylim([-25, 5])
xlabel('Voltage (V)')
xlim([0, 1.2])
%legend()
%title(legend, 'V_{bias} (V)')

%% Do SaP analysis

%Calculate dJ/dV curves
Voc = zeros(num_samples, length(Vbias));
dJdV_Voc = zeros(num_samples, length(Vbias));
h = Vbias(2) - Vbias(1);
for k = 1:num_samples
    for i = 1:length(Vbias)
        J_temp = J_fixed_ion{k,i};
        Voc(k,i) = interp1(J_temp(J_temp~=0), V_fixed_ion(J_temp~=0), 0);
        dJdV = gradient(J_temp(J_temp~=0), V_fixed_ion(J_temp~=0));
        dJdV_Voc(k,i) = interp1(V_fixed_ion(J_temp~=0), dJdV, Voc(k,i));
    end
end

%Extract Vflat

%% Plot dJdV curves
% legend_title = 'N_{ion} (cm^{-3})';
legend_title = 'V_{bi} (V)';
figure('Name', 'SaP-Analysis')

% subplot(1,1,1)
hold on
box on
xline(0, 'black', 'HandleVisibility', 'off')
yline(0, 'black', 'HandleVisibility', 'off')
colours = {[0.9290 0.6940 0.1250], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0 0.4470 0.7410], [0.4940 0.1840 0.5560]};
for k = 1:num_samples
    plot(Vbias, dJdV_Voc(k,:),'Color', colours{k}, 'DisplayName', num2str(1-2*Vbi_offset(k)))%, '%.0e'))
end
xlabel('V_{bias} (V)', 'FontSize', 20)
xlim([Vbias(1), Vbias(end)])
ylabel('dJ/dV|_{Voc} (\Omega^{-1} cm^{-2})', 'FontSize', 20)
legend('Location', 'southeast')
title(legend, legend_title, 'FontSize', 15)


