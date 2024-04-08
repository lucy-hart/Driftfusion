%%
%Plot the flat ion potential versus Vbi for a parameter sweep
%Assumes you have run Ion_Efficiency_v5 before using this

%There's a -1 here as I want to ignore the case with no ions
num_ion_concs = length(solCV(:,1))-1;
num_offsets = length(solCV(1,:));

Symmetric = 1;

Vbi = zeros(1, num_offsets);
TL_offset = zeros(1, num_offsets);
if Symmetric == 1
    Vflat = zeros(num_ion_concs, num_offsets);
elseif Symmetric == 0
    V_ETL = zeros(num_ion_concs, num_offsets);
    V_HTL = zeros(num_ion_concs, num_offsets);
end
V_Emin = zeros(1, num_offsets);
Emin = zeros(1, num_offsets);

for i = 1:num_ion_concs
    for j = 1:num_offsets
        par = solCV{i,j}.par;
        if i == 1
            [V_Emin(j), Emin(j)] = findEminbulk(solCV{end,j});
            Vbi(1,j) = par.Phi_right - par.Phi_left;
            TL_offset(1,j) = par.EF0(end) - par.EF0(1);
        end
        if Symmetric == 1
            Vflat(i,j) = findVflation(solCV{i,j});
        elseif Symmetric == 0
            V_ETL(i,j) = findVinvert_SAM(solCV{i,j}).ETL;
            V_HTL(i,j) = findVinvert_SAM(solCV{i,j}).HTL;
        end
    end
end

%%
figure('Name', 'Vflat and Vbi vs DeltaE', 'Position', [50 50 800 800])
Colours = parula(num_ion_concs);
box on 
for i = 1:num_ion_concs
    hold on
    if i == 1
        plot(Delta_TL, Vbi, 'marker', 'x', 'Color', 'black')
        %plot(Delta_TL, V_Emin, 'marker', 'x', 'Color', 'magenta')
        plot(Delta_TL, TL_offset, 'marker', 'x', 'Color', 'green')
        if Symmetric == 0
            plot(Delta_TL, V_ETL(4,:), 'marker', 'o', 'Color', 'blue', 'LineStyle', 'None', 'MarkerSize', 10)
            plot(Delta_TL, V_HTL(4,:), 'marker', 's', 'Color', 'red', 'LineStyle', 'None', 'MarkerSize', 10)
        end
    end
    if Symmetric == 1
        plot(Delta_TL, Vflat(i,:), 'marker', 'o', 'Color', Colours(i,:))
    end
    plot(Delta_TL, Stats_array(i,:,2), 'marker', 'x', 'Color', Colours(i,:), 'HandleVisibility', 'off')
end
plot(Delta_TL, Stats_array(end,:,2), 'marker', 'x', 'Color', 'black', 'HandleVisibility', 'off')
set(gca, 'Fontsize', 25)
xlabel('Transport Layer Energetic Offset (eV)', 'FontSize', 30)
ylabel('Voltage (V)', 'FontSize', 30)
xlim([0, 0.3])
ylim([0.65, 1.35])
% xticks([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])
% xticklabels({'0.00', '0.05', '0.10', '0.15', '0.20', '0.25', '0.30'})
xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
xticklabels({'0.00', '0.10', '0.20', '0.30', '0.40', '0.50'})
%legend({'V_{BI}', 'V(E_{min,bulk})', '\DeltaE_{F,TL}', 'V_{invert,ETL}', 'V_{invert,HTL}'}, 'Location', 'southwest', 'FontSize', 25)

%%
Run = 0;
    if Run == 1
    figure('Name', 'Vflat and Vbi vs DeltaE', 'Position', [50 50 1600 800])
    Colours = parula(n_ion_concs-1);
    offset = 7;
    box on 
    for i = 1:n_ion_concs
        subplot(1,2,1)
        hold on
        box on
        if i == n_ion_concs
            Vapp = dfana.calcVapp(solCV{i,offset});
            [sample_V, sample_arg] = min(abs(Vapp-V_Emin(offset)));
            x = solCV{i,offset}.x;
            plot(1e7*x, solCV{i,offset}.u(sample_arg-1,:,1), 'Color', 'black', 'HandleVisibility', 'off')
        else
            Vapp = dfana.calcVapp(solCV{i,offset});
            [sample_V, sample_arg] = min(abs(Vapp-Vflat(i,offset)));
            x = solCV{i,offset}.x;
            plot(1e7*x, solCV{i,offset}.u(sample_arg-3,:,1), 'Color', Colours(i,:), 'HandleVisibility', 'off')
        end
    end
    set(gca, 'Fontsize', 25)
    ylabel('Electrostatic Potential (V)', 'FontSize', 30)
    xlabel('Position (nm)', 'FontSize', 30)
    hold off
    
    for i = 1:n_ion_concs-1
        subplot(1,2,2)
        hold on 
        box on
            Vapp = dfana.calcVapp(solCV{i,offset});
            [sample_V, sample_arg] = min(abs(Vapp-Vflat(i,offset)));
            x = solCV{i,offset}.x;
            plot(1e7*x, solCV{i,offset}.u(sample_arg-3,:,4)-solCV{i,offset}.par.Ncat(3), 'Color', Colours(i,:), 'HandleVisibility', 'off')
    end 
    set(gca, 'Fontsize', 25)
    ylabel('Excess Ion Concentration (cm^{-3})', 'FontSize', 30)
    xlabel('Position (nm)', 'FontSize', 30)
    hold off
    end
% xlim([0, 0.3])
% ylim([0.6, 1.45])
% xticks([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])
% xticklabels({'0.00', '0.05', '0.10', '0.15', '0.20', '0.25', '0.30'})
% legend({'V_{BI}', 'V(E_{min,bulk})', 'V_{flat ion}'}, 'Location', 'northeast', 'FontSize', 25)