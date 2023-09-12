%%
%Plot the flat ion potential versus Vbi for a parameter sweep
%Assumes you have run Ion_Efficiency_v5 before using this

%There's a -1 here as I want to ignore the case with no ions
num_ion_concs = length(solCV(:,1))-1;
num_offsets = length(solCV(1,:));

Vbi = zeros(1, num_offsets);
TL_offset = zeros(1, num_offsets);
Vflat = zeros(num_ion_concs, num_offsets);

for i = 1:num_ion_concs
    for j = 1:num_offsets
        par = solCV{i,j}.par;
        if i == 1
            Vbi(1,j) = par.Phi_right - par.Phi_left;
            TL_offset(1,j) = par.EF0(end) - par.EF0(1);
        end
        Vflat(i,j) = findVflation(solCV{i,j});
    end
end

%%
figure('Name', 'Vflat and Vbi vs DeltaE', 'Position', [50 50 800 800])
Colours = parula(n_ion_concs-1);
box on 
for i = 1:num_ion_concs
    hold on
    if i == 1
        plot(Delta_TL, Vbi, 'marker', 'x', 'Color', 'black')
        plot(Delta_TL, TL_offset, 'marker', 'x', 'Color', 'red')
    end
    plot(Delta_TL, Vflat(i,:), 'marker', 'o', 'Color', Colours(i,:))
    plot(Delta_TL, Stats_array(i,:,2), 'marker', 'x', 'Color', Colours(i,:), 'HandleVisibility', 'off')
end
plot(Delta_TL, Stats_array(8,:,2), 'marker', 'x', 'Color', 'black', 'HandleVisibility', 'off')
set(gca, 'Fontsize', 25)
xlabel('Transport Layer Energetic Offset (eV)', 'FontSize', 30)
ylabel('Voltage (V)', 'FontSize', 30)
xlim([0, 0.3])
xticks([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])
xticklabels({'0.00', '0.05', '0.10', '0.15', '0.20', '0.25', '0.30'})
legend({'V_{BI}', '\DeltaE_{F,TLs}', 'V_{flat ion}'}, 'Location', 'northeast', 'FontSize', 25)