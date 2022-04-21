% See how distribution of carriers in PCBM device varies with intensity

%% Read in data
par = pc('Input_files/PTAA_MAPI_PCBM_v2.csv');

%% Do JV scans
eqm = equilibrate(par);
CV_sol_ion_HL = doCV(eqm.ion, 1.0, -0.3, 1.3, -0.3, 1e-3, 1, 321);
CV_sol_ion_LL = doCV(eqm.ion, 0.5, -0.3, 1.3, -0.3, 1e-3, 1, 321);

%% Calculate n values
num_start = sum(CV_sol_ion_HL.par.layer_points(1:2))+1;
num_stop = num_start + CV_sol_ion_HL.par.layer_points(3)-1;
x = CV_sol_ion_HL.par.x_sub;
d = CV_sol_ion_HL.par.d(3);

n_HL = trapz(x(num_start:num_stop), CV_sol_ion_HL.u(:,num_start:num_stop,2),2)/d;
p_HL = trapz(x(num_start:num_stop), CV_sol_ion_HL.u(:,num_start:num_stop,3),2)/d;
n_bar_HL = (n_HL.*p_HL).^0.5;

n_LL = trapz(x(num_start:num_stop), CV_sol_ion_LL.u(:,num_start:num_stop,2),2)/d;
p_LL = trapz(x(num_start:num_stop), CV_sol_ion_LL.u(:,num_start:num_stop,3),2)/d;
n_bar_LL = (n_LL.*p_LL).^0.5;

%% Plotting
x_values = CV_sol_ion_HL.x * 1e7;
colors = {[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],[0 0.4470 0.7410]};

figure(17)
hold on

% Do background colours
% Arguments of patch give the coordinates of the corners of the polygon to
% be shaded
patch('XData',[0 x_values(500) x_values(500) 0], 'YData',[20 20 8 8], 'FaceColor', colors{1},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')
patch('XData',[x_values(500) x_values(600) x_values(600) x_values(500)], 'YData',[20 20 8 8], 'FaceColor', colors{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1100) x_values(1200) x_values(1200) x_values(1100)], 'YData',[20 20 8 8], 'FaceColor', colors{2},...
    'FaceAlpha', 0.5, 'EdgeColor', 'none')
patch('XData',[x_values(1200) x_values(end) x_values(end) x_values(1200)], 'YData',[20 20 8 8],  'FaceColor', colors{3},...
    'FaceAlpha', 0.3, 'EdgeColor', 'none')

plot(x_values, log10(CV_sol_ion_HL.u(138,:,2)), 'b',...
    x_values, log10(CV_sol_ion_HL.u(138,:,3)), 'r',...
    x_values, log10(CV_sol_ion_LL.u(135,:,2)), 'b--', ...
    x_values, log10(CV_sol_ion_LL.u(135,:,3)), 'r--')

hold off

xlim([0, max(x_values)])
ylim([8,20])
ylabel('Charge carrier concentration (cm^{-3})')
xlabel('Device depth (nm)')
legend('','','','','electrons', 'holes', '', '', 'Location', 'bestoutside')
