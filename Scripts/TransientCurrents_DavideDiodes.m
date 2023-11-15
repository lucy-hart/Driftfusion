%File to simulate transients for perovskite/BHJ quasi-tandem devices 
%%
%Set up dark eqm solutions 
read_in_par = 1;

if read_in_par == 1
    parC60 = pc('Input_files/SAM_MAFACsPbIBr_C60.csv');
    parPM6 = pc('Input_files/SAM_MAFACsPbIBr_PM6Y6.csv');
    parPBDBT = pc('Input_files/SAM_MAFACsPbIBr_PBDBTY6.csv');
end

vis = 0;
if vis == 1
    wavelength = 405;
elseif vis == 0
    wavelength = 800;
end

if wavelength == 800
    devices= {parPM6, parPBDBT};
    power = 5;
elseif wavelength == 405
    devices= {parC60, parPM6, parPBDBT};
    power = 10;
end

num_devices = length(devices);
eqm = cell(1, num_devices);

for i = 1:num_devices
   devices{i}.light_source2 = 'laser';
   devices{i}.laser_lambda2 = wavelength;
   devices{i}.pulsepow = power;
   devices{i} = refresh_device(devices{i});
   eqm{i} = equilibrate(devices{i});
end

%%
bias_solution = cell(1, num_devices);
pulse_solution = cell(1, num_devices);
Jt = cell(1, num_devices);

for i = 1:num_devices
    bias_solution{i} = genVappStructs(eqm{i}.ion, -0.5, 1);
    bias_solution{i}.par.int1 = 0;
    bias_solution{i}.par.int2 = 1;
    bias_solution{i}.par.Rs = 10*0.045;
    pulse_solution{i} = doLightPulse(bias_solution{i}, power, 200e-6, 1000, 50, 1, 1);
    Jt{i} = dfana.calcJ(pulse_solution{i}).tot;
end

%%
%Plot 
figure('Name', 'Transient Response')
if wavelength == 405
    colours = {[0 0 0], [1 0 0], [0.4660 0.6740 0.1880]};
elseif wavelength == 800
    colours = {[1 0 0], [0.4660 0.6740 0.1880]};
end

[~,argnorm] = min(abs(pulse_solution{i}.t-100e-6));

box on
hold on
xline(11.3)
yline(0.9)
plot(1e6*pulse_solution{1}.t, ...
    smoothed_square([0, 10, 200e-6, 50],pulse_solution{1}.t)/max(abs(smoothed_square([0, 10, 200e-6, 50],pulse_solution{1}.t))), ...
    'LineStyle', '-', 'HandleVisibility', 'off')
for i = 1:num_devices
    plot(1e6*pulse_solution{i}.t, abs(Jt{i}(:,1)/abs(Jt{i}(argnorm,1))), 'Color', colours{i})
end
hold off 

set(gca, 'FontSize', 25)
xlabel('Time (µs)', 'FontSize', 30)
xlim([0,200])
ylabel('Normalised Current (a.u.)', 'FontSize', 30)
ylim([-0.1, 1.5])
if wavelength == 405
    legend({' C_{60}', ' PM6:Y6', ' PCE12:Y6'}, 'FontSize', 25, 'Location', 'northeast')
elseif wavelength == 800
    legend({' PM6:Y6', ' PCE12:Y6'}, 'FontSize', 25, 'Location', 'northeast')
end