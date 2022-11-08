%SS PL of bi-layers with different thicknesses of C60
%Illuminate through the C60

par1 = pc('Input_files/MAPbI_C60_30nm.csv');
par2 = pc('Input_files/MAPbI_C60_60nm.csv');
par3 = pc('Input_files/MAPbI_C60_90nm.csv');

pars = {par1, par2, par3};
num_films = length(pars);
eqm = cell(1,num_films);

for i = 1:num_films
    par = pars{i};
    par.laser_lambda2 = 532;
    par.pulsepow = 50;
    pars{i} = refresh_device(par);
    eqm{i} = equilibrate(pars{i});
end
%% Find solutoions at different intensisites 
num_samples = 1;
illuminated_sols = cell(num_samples, num_films);
%LightInt = logspace(log10(0.1), log10(5), num_samples-1);
%LightInt = [LightInt, 1];
LightInt = [1];

for i = 1:num_films
    for j = 1:num_samples
        illuminated_sols{j,i} = changeLight(eqm{i}.ion, LightInt(j), 0, 2);
    end
end

%% Calculate Rad, Bulk and Surface recombination
Rad = zeros(num_samples, num_films);
Bulk = zeros(num_samples, num_films);
Surface = zeros(num_samples, num_films);

for i = 1:num_films
    for j = 1:num_samples
        sol = illuminated_sols{j,i};
        loss_currents = dfana.calcr(sol,'sub');
        x = sol.par.x_sub;
        num_points = length(x);
        J = dfana.calcJ(sol);

        Rad(j, i) = trapz(x, loss_currents.btb(end,:), 2);
        Bulk(j,i) = trapz(x, loss_currents.srh(end,:), 2);
        Surface(j,i) = trapz(x, loss_currents.vsr(end,:), 2);
    end
end
%%
x = [1,2,3];
num = 1;
y = [Rad(num,1) Rad(num,2) Rad(num,3); ...
    Bulk(num,1) Bulk(num,2) Bulk(num,3);...
    Surface(num,1), Surface(num,2), Surface(num,3)];
markers = {'x', 'o', 's'};

figure('Name', 'Recombination','Position', [100 100 1250 2000])
box on

for i = 1:num_films
    semilogy(x, y(:,i), 'LineStyle', 'none', 'Marker', markers{i}, 'MarkerSize', 10)
    hold on
end

set(gca, 'FontSize', 30)
xlim([0, 4])
xticks([1,2,3])
xticklabels({'Radiative', 'Bulk', 'Surface'})
ylabel('Recombination Rate (s^{-1})', 'FontSize', 30)
legend({'30 nm', '60 nm', '90 nm'}, 'Location', 'northeast', 'FontSize', 30, 'Location', 'northwest')
title(legend, 'C_{60} Thickness')
%%
figure('Name', 'PL vs LI')
B = eqm{1}.ion.par.B(3);
for i = 1:num_films
    loglog(LightInt, B*Rad(:,i))
    hold on 
end
xlabel('Intensity')
xlim([0.1, 5])
ylabel('Radiative Recombination (s^{-1})')
legend({'30 nm', '60 nm', '90 nm'}, 'Location', 'northwest')
title(legend, 'C_{60} Thickness')