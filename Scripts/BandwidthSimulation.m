%File to simulate bandwidth measurements
%%
%Set up dark eqm solutions 
read_in_par = 1;

if read_in_par == 1
    parC60 = pc('Input_files/SAM_MAFACsPbIBr_C60_EunyoungValues.csv');
end

wavelength = 532;
devices= {parC60};
power = 2;
frequency = [10 50 100 500 1000 5000];
num_cycles = 30;

QandD = 0;

num_devices = length(devices);
num_f = length(frequency);
eqm = cell(1, num_devices);

for i = 1:num_devices
   devices{i}.light_source2 = 'laser';
   devices{i}.laser_lambda2 = wavelength;
   devices{i}.pulsepow = power;
   devices{i} = refresh_device(devices{i});
   eqm{i} = equilibrate(devices{i});
end

%%
bias_solution = cell(num_f, num_devices);
pulse_solution = cell(num_f, num_devices);
Jt = cell(num_f, num_devices);

for i = 1:num_devices
    for j = 1:num_f
        period = 1/frequency(j);
        if QandD == 0
            bias_solution{j,i} = genVappStructs(eqm{i}.ion, -0.5, 1);
        else
            bias_solution{j,i} = eqm{i}.ion;
        end
        bias_solution{j,i}.par.int1 = 0;
        bias_solution{j,i}.par.int2 = 1;
        %bias_solution{j,i}.par.Rs = 20*0.045;
        pulse_solution{j,i} = doSinWave(bias_solution{j,i}, power, num_cycles*period, 5000, frequency(j), 1);
        Jt{j,i} = dfana.calcJ(pulse_solution{j,i}).tot;
    end
end

%%
%Plot 
figure('Name', 'Transient Response')
colours = {[0 0 0], [1 0 0], [0.4660 0.6740 0.1880]};

box on
hold on
%xline(11.3)
%yline(0.9)
for i = 1:num_f
    plot(pulse_solution{i,1}.t, Jt{i,1}(:,1), 'Color', colours{i})
end
hold off 

set(gca, 'FontSize', 25)
xlabel('Time (s)', 'FontSize', 30)
%xlim([0,200])
ylabel('Current (a.u.)', 'FontSize', 30)
%ylim([-0.1, 1.5])
