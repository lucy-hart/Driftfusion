%par = pc('Input_files/PTAA_MAPI_NoOffset_FastJVPaperParams.csv');
par = pc('Input_files/PTAA_MAPI_NoOffset.csv');
par.RelTol_vsr = 0.075;

%IonConc = logspace(16,18,3);
IonConc = [5e16 1e17 5e17];
%IonConc = 1e19;
num_samples = length(IonConc);
taus = [20e-9 100e-9 200e-9];
%mobs = flip(logspace(-1.3,0,num_samples));
devices = cell(num_samples, 1);

results_ion = zeros(num_samples, 4);
results_bias = zeros(num_samples, 4);
results_el = zeros(1,4);

j_ion = cell(1,num_samples);
j_bias = cell(1, num_samples);

for i = 1:num_samples
    par.Ncat(:) = IonConc(i);
    par.Nani(:) = IonConc(i);
%     par.mu_n(3) = mobs(i);
%     par.mu_p(3) = mobs(i);
    par.taun(3) = taus(i);
    par.taup(3) = taus(i);
    par = refresh_device(par);
    devices{i} = equilibrate(par);
end

%%
suns = 1;

for i = 1:num_samples
    try
        JV_sol_ion = doCV(devices{i}.ion, suns, -0.2, 1.2, -0.2, 1e-4, 1, 281);
    catch
        try
            JV_sol_ion = doCV(devices{i}.ion, suns, -0.2, 1.18, -0.2, 1e-4, 1, 277);
            Vapp2 = dfana.calcVapp(JV_sol_ion);
        catch
            JV_sol_ion = doCV(devices{i}.ion, suns, -0.2, 1.15, -0.2, 1e-4, 1, 271);
            Vapp3 = dfana.calcVapp(JV_sol_ion);
        end
    end
    
    stats_ion = CVstats(JV_sol_ion);
    results_ion(i,1) = -stats_ion.Jsc_f;
    results_ion(i,2) = stats_ion.Voc_f;
    results_ion(i,3) = stats_ion.FF_f;
    results_ion(i,4) = stats_ion.efficiency_f;
           
    j_ion{i} = dfana.calcJ(JV_sol_ion).tot(:,1);

    %Fix this value based on the Voc of the default parameters (Nion = 5e16
    %cm-3) - better reflects the paper where they don't change scan range
    %to take account Voc decreases during aging
    V_bias = 1.15;
    
    biased_eqm_ion = genVappStructs(devices{i}.ion, V_bias, 1);
    illuminated_sol_ion = changeLight(biased_eqm_ion, suns, 0, 1);
    illuminated_sol_ion.par.mobseti = 0;
    
    try
        JV_sol_bias = doCV(illuminated_sol_ion, suns, -0.2, 1.2, -0.2, 0.1, 1, 281);
    catch
        try
            JV_sol_bias = doCV(illuminated_sol_ion, suns, -0.2, 1.18, -0.2, 0.1, 1, 277);
            Vapp2 = dfana.calcVapp(JV_sol_ion);
        catch
            JV_sol_bias = doCV(illuminated_sol_ion, suns, -0.2, 1.15, -0.2, 0.1, 1, 271);
            Vapp3 = dfana.calcVapp(JV_sol_ion);
        end
    end
    
    stats_bias = CVstats(JV_sol_bias);
    results_bias(i,1) = -stats_bias.Jsc_f;
    results_bias(i,2) = stats_bias.Voc_f;
    results_bias(i,3) = stats_bias.FF_f;
    results_bias(i,4) = stats_bias.efficiency_f;

    j_bias{i} = dfana.calcJ(JV_sol_bias).tot(:,1);

end

JV_sol_el = doCV(devices{1}.el, suns, -0.2, 1.2, -0.2, 0.1, 1, 281);

stats_el = CVstats(JV_sol_el);
results_el(1) = -stats_el.Jsc_f;
results_el(2) = stats_el.Voc_f;
results_el(3) = stats_el.FF_f;
results_el(4) = stats_el.efficiency_f;

j_el = dfana.calcJ(JV_sol_el).tot(:,1);

Vapp = dfana.calcVapp(JV_sol_el);

%%
figure('Name', 'StatsvsNion')

for i = 1:num_samples
    for j = 1:4
        subplot(2,2,j)
        hold on
        if i == 1
            plot(0, abs(results_el(j)), 'Marker', 'o', 'Color', 'k')
        end
        plot(i, abs(results_ion(i,j)), 'Marker', 'square', 'Color', 'blue')
        plot(i, abs(results_bias(i,j)), 'Marker', 'diamond', 'Color', 'red')
        hold off
    end
end

%%
figure('Name', 'JVs', 'Position', [100 100 1500 750])
cmap_ions = colormap(autumn(num_samples+1));
cmap_bias = colormap(winter(num_samples+1));
labels = {'1.0 x 10^{16}', '3.2 x 10^{16}', '1.0 x 10^{17}', '3.2 x 10^{17}', '1.0 x 10^{18}', 'None'};

for i = 1:num_samples
        hold on
        subplot(1,2,1)
        try
            plot(Vapp, 1e3*j_ion{i}, 'Color', cmap_ions(i+1,:))
        catch 
            num_points = length(j_ion{i});
            if length(Vapp2) == num_points
                plot(Vapp2, 1e3*j_ion{i}, 'Color', cmap_ions(i+1,:))
            else
                plot(Vapp3, 1e3*j_ion{i}, 'Color', cmap_ions(i+1,:))
            end
        end
        hold off

        hold on
        subplot(1,2,2)
        try
            plot(Vapp, 1e3*j_bias{i}, 'Color', cmap_bias(i+1,:))
        catch 
            num_points = length(j_bias{i});
            if length(Vapp2) == num_points
                plot(Vapp2, 1e3*j_bias{i}, 'Color', cmap_bias(i+1,:))
            else
                plot(Vapp3, 1e3*j_bias{i}, 'Color', cmap_bias(i+1,:))
            end
        end
        hold off
end

subplot(1,2,1)
hold on
plot(Vapp, 1e3*j_el, 'Color', 'k')
xline(0, 'black', 'LineWidth', 1)
yline(0, 'black', 'LineWidth', 1)
hold off
set(gca, 'FontSize', 25)
ylim([-23, 3])
xlim([-0.1, 1.2])
xlabel('Voltage (V)', 'FontSize', 30)
ylabel('Current Density (mA cm^{-2})', 'FontSize', 30)
legend(labels, 'FontSize', 20, 'Location', 'northwest', 'NumColumns', 2)
title(legend, 'Mobile Ion Density (cm^{-3})', 'FontSize', 20)

subplot(1,2,2)
hold on
plot(Vapp, 1e3*j_el, 'Color', 'k')
xline(0, 'black', 'LineWidth', 1)
yline(0, 'black', 'LineWidth', 1)
hold off
set(gca, 'FontSize', 25)
ylim([-23, 3])
xlim([-0.1, 1.2])
xlabel('Voltage (V)', 'FontSize', 30)
ylabel('Current Density (mA cm^{-2})', 'FontSize', 30)
legend(labels, 'FontSize', 20, 'Location', 'northwest', 'NumColumns', 2)
title(legend, 'Mobile Ion Density (cm^{-3})', 'FontSize', 20)
