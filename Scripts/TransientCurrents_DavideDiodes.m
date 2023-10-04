%File to simulate transients for perovskite/BHJ quasi-tandem devices 
%%
%Set up dark eqm solutions 
parC60 = pc('Input_files/SAM_MAFACsPbIBr_C60.csv');
parPM6 = pc('Input_files/SAM_MAFACsPbIBr_PM6Y6.csv');

devices= {parC60, parPM6};
num_devices = length(devices);
eqm = cell(1, num_devices);

wavelength = 405;

for i = 1:num_devices
   devices{i}.light_source2 = 'laser';
   devices{i}.laser_lambda2 = wavelength;
   eqm{1,i} = equilibrate(devices{i});
end

%%
bias_solution = cell(1, num_devices);
pulse_solution = cell(1, num_devices);
Jt = cell(1, num_devices);
for i = 1:num_devices
    bias_solution{i} = genVappStructs(eqm{i}.ion, -0.5, 1);
    pulse_solution{i} = doLightPulse(SC_solution{i},10, 200e-6, 1000, 50, 1, 0);
    Jt{i} = dfana.calcJ(pulse_solution{i}).tot;
end

%%
%Plot 
figure('Name', 'Linear Dynamic Response')

box on
hold on
for i = 1:num_devices
    plot(1e6*pulse_solution{i}.t, Jt{i}(:,1)/max(abs(Jt{i}(:,1))))
end
hold off 

set(gca, 'FontSize', 25)
xlabel('Time (\mus)', 'FontSize', 30)
ylabel('Normalised Current (a.u.)', 'FontSize', 30)

legend({'C_{60}', ' PM6:Y6', ' PM7:Y6', ' PCE12:Y6'}, 'FontSize', 25, 'Location', 'northwest')