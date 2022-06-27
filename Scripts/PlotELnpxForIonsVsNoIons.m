% Plot energy levels and carrier popultions at short circuit 
% with and without ion motion.
% SC conditions for forward scan are at 31st timepoint, and 291st for reverse 
% OC at 106, 137 and 183 for Kloc-6, PCBM and ICBA respectively

num = 2;
%T = 31;
T = OC_time_ion(num);
lower_lim=[-7.5, -7.5, -7.5];

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

T = OC_time_el(num);

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
 semilogy(x_values, 18*ones(1,length(x_values)), 'black','LineStyle', '--')

hold off

xlim([0, max(x_values)])
ylim([9,19])
ylabel('Charge carrier concentration (cm^{-3})')
xlabel('Device depth (nm)')
legend('','','','','electrons', 'holes', '', '', 'cations','','', 'Location', 'bestoutside')

%% Plot el vs ion populations at OC vs SC

num = 4;
T = 31;
x_values = CV_solutions_ion{num}.x * 1e7;

figure('Name', 'Carrier Densities', 'Position', [100 100 1000 2000])

subplot(3,1,1)
hold on

% Do background colours
% Arguments of patch give the coordinates of the corners of the polygon to
% be shaded
patch('XData',[0 x_values(500) x_values(500) 0], 'YData',[30 30 2 2], 'FaceColor', colors_JV{1},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
patch('XData',[x_values(500) x_values(600) x_values(600) x_values(500)], 'YData',[30 30 2 2], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1100) x_values(1200) x_values(1200) x_values(1100)], 'YData',[30 30 2 2], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1200) x_values(end) x_values(end) x_values(1200)], 'YData',[30 30 2 2],  'FaceColor', colors_JV{3},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')

semilogy(x_values, log10(CV_solutions_ion{num}.u(T,:,2)), 'b-',...
    x_values, log10(CV_solutions_ion{num}.u(T,:,3)), 'r-',...    
    x_values, log10(CV_solutions_el{num}.u(T,:,2)), 'b--',...
    x_values, log10(CV_solutions_el{num}.u(T,:,3)), 'r--')

ylim([9,20])
ylabel('Carrier Density (cm^{-3})')
xlim([0, max(x_values)])
xlabel('Device depth (nm)')
legend('','','','','electrons','holes', 'NumColumns', 1, 'Position', [0.7,0.85,0.05,0.03])

hold off 

subplot(3,1,2)
T = OC_time_ion(num);

hold on

% Do background colours
% Arguments of patch give the coordinates of the corners of the polygon to
% be shaded
patch('XData',[0 x_values(500) x_values(500) 0], 'YData',[30 30 2 2], 'FaceColor', colors_JV{1},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
patch('XData',[x_values(500) x_values(600) x_values(600) x_values(500)], 'YData',[30 30 2 2], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1100) x_values(1200) x_values(1200) x_values(1100)], 'YData',[30 30 2 2], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1200) x_values(end) x_values(end) x_values(1200)], 'YData',[30 30 2 2],  'FaceColor', colors_JV{3},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')

semilogy(x_values, log10(CV_solutions_ion{num}.u(T,:,2)), 'b-',...
         x_values, log10(CV_solutions_ion{num}.u(T,:,3)), 'r-')

T = OC_time_el(num);

semilogy(x_values, log10(CV_solutions_el{num}.u(T,:,2)), 'b--',...
         x_values, log10(CV_solutions_el{num}.u(T,:,3)), 'r--')

hold off 

xlim([0, max(x_values)])
ylim([9,20])
ylabel('Carrier Density (cm^{-3})')
xlabel('Device depth (nm)')

subplot(3,1,3)
T = 31;
upper_lim_ion = [4e18, 4e18, 3.5e18, 3e18];

hold on

% Do background colours
% Arguments of patch give the coordinates of the corners of the polygon to
% be shaded
patch('XData',[0 x_values(500) x_values(500) 0], 'YData',[1e19 1e19 1e9 1e9], 'FaceColor', colors_JV{1},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
patch('XData',[x_values(500) x_values(600) x_values(600) x_values(500)], 'YData',[1e19 1e19 1e9 1e9], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1100) x_values(1200) x_values(1200) x_values(1100)], 'YData',[1e19 1e19 1e9 1e9], 'FaceColor', colors_JV{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1200) x_values(end) x_values(end) x_values(1200)], 'YData',[1e19 1e19 1e9 1e9],  'FaceColor', colors_JV{3},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')

plot(x_values, (CV_solutions_ion{num}.u(T,:,4)), 'color', colors_JV{2})

T = OC_time_ion(num);

 semilogy(x_values, (CV_solutions_ion{num}.u(T,:,4)), 'color', colors_JV{2}, 'LineStyle', '--')

hold off 

xlim([0, max(x_values)])
ylim([1e10, upper_lim_ion(num)])
ylabel('Cation Population (cm^{-3})')
xlabel('Device depth (nm)')
legend('','','','','Short Circuit','Open Circuit', 'NumColumns', 1, 'Position', [0.7,0.25,0.05,0.03])

fig = gcf;

%%
%exportgraphics(fig, ...
 %   'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Weidong_ETL\Simulations\v2\Carrier_Distributions_PCBM_el_and_ion_BIG.png', ...
  %  'Resolution', 300)
