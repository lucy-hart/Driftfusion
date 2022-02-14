% Plot energy levels and carrier popultions at short circuit 
% with and without ion motion.
% SC conditions for forward scan are at 31st timepoint, and 291st for reverse 
% OC at 106, 137 and 183 for Kloc-6, PCBM and ICBA respectively

num = 1;
T = 31;
%T = OC_time_ion(num);
lower_lim=[-6.5, -7, -7.5];

[Ecb_ion, Evb_ion, Efn_ion, Efp_ion] = dfana.calcEnergies(CV_solutions_ion{num});
[Ecb_el, Evb_el, Efn_el, Efp_el] = dfana.calcEnergies(CV_solutions_el{num});
x_values = CV_solutions_ion{num}.x * 1e7;

figure(5)

subplot(2,1,1)
hold on

% Do background colours
% Arguments of patch give the coordinates of the corners of the polygon to
% be shaded
patch('XData',[0 x_values(500) x_values(500) 0], 'YData',[-2 -2 lower_lim(num) lower_lim(num)], 'FaceColor', colors_JV{1},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
patch('XData',[x_values(500) x_values(600) x_values(600) x_values(500)], 'YData',[-2 -2 lower_lim(num) lower_lim(num)], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1100) x_values(1200) x_values(1200) x_values(1100)], 'YData',[-2 -2 lower_lim(num) lower_lim(num)], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1200) x_values(end) x_values(end) x_values(1200)], 'YData',[-2 -2 lower_lim(num) lower_lim(num)],  'FaceColor', colors_JV{3},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')

plot(x_values, Ecb_ion(T,:), 'b',...
    x_values, Evb_ion(T,:), 'r',...
    x_values, Efn_ion(T,:), 'b--', ...
    x_values, Efp_ion(T,:), 'r--')

hold off

xlim([0, max(x_values)])
ylim([lower_lim(num), -2.0])
ylabel('Energy, eV')
xlabel('Device depth (nm)')
%title('Energy Levels with Ions')

subplot(2,1,2)

hold on

% Do background colours
% Arguments of patch give the coordinates of the corners of the polygon to
% be shaded
patch('XData',[0 x_values(500) x_values(500) 0], 'YData',[-2 -2 lower_lim(num) lower_lim(num)], 'FaceColor', colors_JV{1},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
patch('XData',[x_values(500) x_values(600) x_values(600) x_values(500)], 'YData',[-2 -2 lower_lim(num) lower_lim(num)], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1100) x_values(1200) x_values(1200) x_values(1100)], 'YData',[-2 -2 lower_lim(num) lower_lim(num)], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1200) x_values(end) x_values(end) x_values(1200)], 'YData',[-2 -2 lower_lim(num) lower_lim(num)],  'FaceColor', colors_JV{3},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')

%T = OC_time_el(num);

plot(x_values, Ecb_el(T,:), 'b',...
    x_values, Evb_el(T,:), 'r',...
    x_values, Efn_el(T,:), 'b--', ...
    x_values, Efp_el(T,:), 'r--')

hold off 

xlim([0, max(x_values)])
ylim([lower_lim(num), -2.0])
ylabel('Energy, eV')
xlabel('Device depth (nm)')
legend('','','','','E_{C}', 'E_{V}', 'E_{f,n}', 'E_{f,p}', 'NumColumns', 2, 'Location', 'northeast')
%title('Energy levels without Ions')

%%
figure(6)

num = 3;
T=31;

hold on

% Do background colours
% Arguments of patch give the coordinates of the corners of the polygon to
% be shaded
patch('XData',[0 x_values(500) x_values(500) 0], 'YData',[19 19 9 9], 'FaceColor', colors_JV{1},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
patch('XData',[x_values(500) x_values(600) x_values(600) x_values(500)], 'YData',[19 19 9 9], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1100) x_values(1200) x_values(1200) x_values(1100)], 'YData',[19 19 9 9], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1200) x_values(end) x_values(end) x_values(1200)], 'YData',[19 19 9 9],  'FaceColor', colors_JV{3},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')

semilogy(x_values, log10(CV_solutions_ion{num}.u(OC_time_ion(num),:,2)), 'b-',...
    x_values, log10(CV_solutions_ion{num}.u(OC_time_ion(num),:,3)), 'r-',...    
    x_values, log10(CV_solutions_ion{num}.u(T,:,2)), 'b--',...
    x_values, log10(CV_solutions_ion{num}.u(T,:,3)), 'r--')

 semilogy(x_values, log10(CV_solutions_ion{num}.u(OC_time_ion(num),:,4)), 'color', colors_JV{2})
 semilogy(x_values, log10(CV_solutions_ion{num}.u(T,:,4)), 'color', colors_JV{2}, 'LineStyle', '--')
 semilogy(x_values, 17*ones(1,length(x_values)), 'black','LineStyle', '--')

hold off

xlim([0, max(x_values)])
ylim([9,19])
ylabel('Charge carrier concentration (cm^{-3})')
xlabel('Device depth (nm)')
legend('','','','','electrons', 'holes', '', '', 'cations','','', 'Location', 'bestoutside')

%
%x_values, log10(CV_solutions_ion{num}.u(T,:,4)), 'color', colors_JV{2}, 'LineStyle', '--', ...
