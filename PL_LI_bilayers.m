% Code to compare light intensity dependence of QFLS in bilayer device
% stacks with ETMs with varying LUMOs
AllSamples = 0;
if AllSamples == 0
    filenames = ["Input_files/FACsPbIBr_PCBM_Bilayer.csv"];
elseif AllSamples == 1
    filenames = ["Input_files/FACsPbIBr_PCBM_Bilayer_lowerLUMO.csv",...
    "Input_files/FACsPbIBr_PCBM_Bilayer.csv",...
    "Input_files/FACsPbIBr_PCBM_Bilayer_higherLUMO.csv"];
end

num_samples = length(filenames);
samples = cell(1, num_samples);
eqm = cell(1, num_samples);
for i = 1:num_samples
    par_bi = pc(filenames(i));
    par_bi.light_source1 = 'laser';
    par_bi.laser_lambda1 = 532; 
    par_bi.pulsepow = 100; %I think this is in mAcm-2 so I have used 100 (same intensity as AM1.5G)
    samples{i} = refresh_device(par_bi);
    eqm{i} = equilibrate(samples{i});
end

%% Find illuminated SS solutions
num_int = 10;
illuminated_sol_bilayer = cell(num_int,num_samples);
intensities = logspace(-1, 1, num_int);

for i = 1:num_int
    for j = 1:num_samples
    illuminated_sol_bilayer{i,j} = lightonRs(eqm{j}.ion, intensities(i), 5e-3, 1, 1e6, 1000);
    end
end 
%% Get PL
PL = zeros(num_int, num_samples);
for i = 1:num_int
    for j = 1:num_samples
    sol = illuminated_sol_bilayer{i,j};
    loss_currents = dfana.calcr(sol,'sub');
    x = sol.par.x_sub;
    num_points = length(x);
    PL(i,j) = trapz(x, loss_currents.btb(end,:), 2)';
    end
end    

%% Plot
colours = {[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],[0 0.4470 0.7410]};
figure('Name', 'PL vs LI')

for i = 1: num_samples
    loglog(intensities, PL(:,i), 'color', colours{i}, 'marker', 'o')
    hold on 
end 
xlabel('Intensity')
ylabel('PL Counts')

figure('Name', 'PL vs LI Gradients')

for i = 1: num_samples
    semilogx(intensities, gradient(log(PL(:,i)), log(intensities)), 'color', colours{i}, 'marker', 'x')
    hold on 
end 
xlabel('Intensity')
ylabel('Gradient log(PL) vs log(Intensity)')
