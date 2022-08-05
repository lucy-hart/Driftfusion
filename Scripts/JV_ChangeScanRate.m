% Make plot of JV curves at different scan rates

num = 1;
scan_rates = [1e-2, 1e-1, 1];
JV_Solutions = cell(1,5);

j=1;

for i = scan_rates
    try
        JV_Solutions{j} = doCV(eqm_solutions_dark{num}.ion, 1.15, -0.2, 1.2, -0.2, i, 1, 301);
    catch
        warning(['JV scan failed for scan rate ' num2str(i) ' Vs-1.'])
    end
    j = j+1;
end

%% Plotting 

figure('Name', 'Scan Rate Dependent JV', 'Position', [100 100 1250 1250])
for k = 1:3
plot(v(1:151), 1000*dfana.calcJ(JV_Solutions{k}).tot(1:151,1), 'color', colors_JV{k})
hold on
plot(v(151:end), 1000*dfana.calcJ(JV_Solutions{k}).tot(151:end,1), 'color', colors_JV{k}, 'LineStyle', '--')
hold on
end

plot(v(:), zeros(1,length(v)), 'black', 'LineWidth', 1)
hold off
legend({'0.01 Vs^{-1}','','0.1 Vs^{-1}','', '1 Vs^{-1}',''}, 'Location', 'northwest')
xlim([-0.1, 1.2])
ylim([-25, 5])
xlabel('Voltage(V)')
ylabel('Current Density (Acm^{-2})')