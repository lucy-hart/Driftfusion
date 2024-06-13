MyParams = 1;
if MyParams == 0
    par = pc('Input_files/PTAA_MAPI_NoOffset_FastJVPaperParams.csv');
    v_surf_l = logspace(log10(500), log10(50000), 5);
    v_surf_r = logspace(log10(2000), log10(200000), 5);
    labels = {'2000', '6300', '20000', '63000', '200000'};
    %Fix this value based on the Voc of the default parameters - better reflects 
    % the paper where they don't change scan range
    %to take account Voc decreases during aging
    V_bias = 1.15;
elseif MyParams == 1
    par = pc('Input_files/PTAA_MAPI_NoOffset.csv');
    v_surf_l = logspace(log10(5), log10(500), 5);
    v_surf_r = logspace(log10(20), log10(2000), 5);
    labels = {'20', '63', '200', '630', '2000'};
    %See comment above
    V_bias = 1.157;
end
par.RelTol_vsr = 0.075;

num_samples = length(v_surf);
devices = cell(num_samples, 1);

results_ion = zeros(num_samples, 4);
results_bias = zeros(num_samples, 4);
results_el = zeros(num_samples,4);

j_ion = cell(1,num_samples);
j_bias = cell(1, num_samples);
j_el = cell(1, num_samples);

for i = 1:num_samples
    par.sn(2) = v_surf_l(i);
    par.sp(4) = v_surf_r(i);
    par = refresh_device(par);
    devices{i} = equilibrate(par);
end

%%
suns = 1;

for i = 1:num_samples
    try
        JV_sol_ion = doCV(devices{i}.ion, suns, -0.2, 1.2, -0.2, 1e-4, 1, 281);
        Vapp = dfana.calcVapp(JV_sol_ion);
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

    try
        JV_sol_el = doCV(devices{i}.el, suns, -0.2, 1.2, -0.2, 1e-4, 1, 281);
        Vapp = dfana.calcVapp(JV_sol_el);
    catch
        try
            JV_sol_el = doCV(devices{i}.el, suns, -0.2, 1.18, -0.2, 1e-4, 1, 277);
            Vapp2 = dfana.calcVapp(JV_sol_el);
        catch
            JV_sol_el = doCV(devices{i}.el, suns, -0.2, 1.15, -0.2, 1e-4, 1, 271);
            Vapp3 = dfana.calcVapp(JV_sol_el);
        end
    end
    
    stats_el = CVstats(JV_sol_el);
    results_el(i,1) = -stats_el.Jsc_f;
    results_el(i,2) = stats_el.Voc_f;
    results_el(i,3) = stats_el.FF_f;
    results_el(i,4) = stats_el.efficiency_f;
           
    j_el{i} = dfana.calcJ(JV_sol_el).tot(:,1);
    
    biased_eqm_ion = genVappStructs(devices{i}.ion, V_bias, 1);
    illuminated_sol_ion = changeLight(biased_eqm_ion, suns, 0, 1);
    illuminated_sol_ion.par.mobseti = 0;
    
    try
        JV_sol_bias = doCV(illuminated_sol_ion, suns, -0.2, 1.2, -0.2, 0.1, 1, 281);
        Vapp = dfana.calcVapp(JV_sol_ion);
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

%%
figure('Name', 'StatsvsNion')

for i = 1:num_samples
    for j = 1:4
        subplot(2,2,j)
        hold on
        plot(i, abs(results_el(i,j)), 'Marker', 'o', 'Color', 'black')
        plot(i, abs(results_ion(i,j)), 'Marker', 'square', 'Color', 'blue')
        plot(i, abs(results_bias(i,j)), 'Marker', 'diamond', 'Color', 'red')
        hold off
    end
end

%%
figure('Name', 'JVs', 'Position', [100 100 2000 750])
cmap_ions = colormap(autumn(num_samples+1));
cmap_bias = colormap(winter(num_samples+1));
cmap_el = colormap(summer(num_samples+1));

for i = 1:num_samples
        hold on
        subplot(1,3,1)
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
        subplot(1,3,2)
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

        hold on
        subplot(1,3,3)
        try
            plot(Vapp, 1e3*j_el{i}, 'Color', cmap_el(i+1,:))
        catch 
            num_points = length(j_el{i});
            if length(Vapp2) == num_points
                plot(Vapp2, 1e3*j_el{i}, 'Color', cmap_el(i+1,:))
            else
                plot(Vapp3, 1e3*j_el{i}, 'Color', cmap_el(i+1,:))
            end
        end
        hold off
end

subplot(1,3,1)
hold on
xline(0, 'black', 'LineWidth', 1)
yline(0, 'black', 'LineWidth', 1)
hold off
set(gca, 'FontSize', 25)
ylim([-23, 3])
xlim([-0.1, 1.2])
xlabel('Voltage (V)', 'FontSize', 30)
ylabel('Current Density (mA cm^{-2})', 'FontSize', 30)
legend(labels, 'FontSize', 20, 'Location', 'northwest')
title(legend, 'v_{surf,ETL} (cm s^{-1})', 'FontSize', 20)

subplot(1,3,2)
hold on
xline(0, 'black', 'LineWidth', 1)
yline(0, 'black', 'LineWidth', 1)
hold off
set(gca, 'FontSize', 25)
ylim([-23, 3])
xlim([-0.1, 1.2])
xlabel('Voltage (V)', 'FontSize', 30)
ylabel('Current Density (mA cm^{-2})', 'FontSize', 30)
legend(labels, 'FontSize', 20, 'Location', 'northwest')
title(legend, 'v_{surf,ETL} (cm s^{-1})', 'FontSize', 20)

subplot(1,3,3)
hold on
xline(0, 'black', 'LineWidth', 1)
yline(0, 'black', 'LineWidth', 1)
hold off
set(gca, 'FontSize', 25)
ylim([-23, 3])
xlim([-0.1, 1.2])
xlabel('Voltage (V)', 'FontSize', 30)
ylabel('Current Density (mA cm^{-2})', 'FontSize', 30)
legend(labels, 'FontSize', 20, 'Location', 'northwest')
title(legend, 'v_{surf,ETL} (cm s^{-1})', 'FontSize', 20)