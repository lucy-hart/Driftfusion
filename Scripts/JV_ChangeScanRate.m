% Make plot of JV curves at different scan rates

num = 2;
scan_rates = [1e-3, 1e-2, 0.05];
JV_Solutions = cell(1,4);

j=1;

for i = scan_rates
    JV_Solutions{j} = doCV(eqm_solutions_dark{num}.ion, 1.15, -0.3, 1.25, -0.3, i, 1, 321);
    j = j+1;
end

%% Plotting 

figure(4)
for k = 1:3
plot(v(1:161), -dfana.calcJ(JV_Solutions{k}).tot(1:161,1), 'color', colors_JV{k})
hold on
plot(v(161:end), -dfana.calcJ(JV_Solutions{k}).tot(161:end,1), 'color', colors_JV{k}, 'LineStyle', '--')
hold on
end

plot(v(:), zeros(1,length(v)), 'black', 'LineWidth', 1)
hold off
legend({'1 mVs^{-1}','','10 mVs^{-1}','', '50 mVs^{-1}',''}, 'Location', 'southwest')
xlim([0, 1.25])
ylim([0, 0.025])
xlabel('Voltage(V)')
ylabel('Current Density (Acm^{-2})')