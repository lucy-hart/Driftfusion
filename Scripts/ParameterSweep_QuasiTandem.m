%Use this file to symmetrically sweep HOMO/LUMO offsets vs ion
%concentration

%TURN SAVE OFF TO START OFF WITH (final cell)

tic
%% Define parameter space
Donor_HOMO = linspace(-5.4, -5.1, 13);
k_rec = logspace(-11, -10, 11);
n_k_rec = length(k_rec); 
n_HOMOs = length(Donor_HOMO);

params = cell(n_k_rec, n_HOMOs);

for i=1:n_k_rec
    for j=1:n_HOMOs
            params{i,j} = [k_rec(i), Donor_HOMO(j)];
    end
end

%%
error_log = zeros(n_k_rec, n_HOMOs);
soleq = cell(n_k_rec, n_HOMOs);
solCV = cell(n_k_rec, n_HOMOs);
results = cell(n_k_rec, n_HOMOs);

%% Do (many) JV sweeps
%Remeber to update the work function values if you change these parameters
%between files 
illumination = 1;

%Reset the electrode work functions in each loop to be safe as they are
%changed for the cases where E_LUMO (E_HOMO) is far below (above) the CB
%(VB)
for i = 1:n_k_rec
    for j = 1:n_HOMOs
        disp(["Donor_HOMO = ", num2str(Donor_HOMO(j)), " eV"])
        disp(["k_rec = ", num2str(k_rec(i)), " cm{^-3}"])
        
        par = pc('Input_files/SAM_MAFACsPbIBr_PM6Y6.csv');

        %Donor HOMO Energetics
        par.Phi_IP(5) = params{i,j}(2);
        par.EF0(5) = (par.Phi_IP(5)+par.Phi_EA(5))/2;
        par.Et(5) = (par.Phi_IP(5)+par.Phi_EA(5))/2;

        %BHJ recombination rate
        par.B(5) = params{i,j}(1);

        par = refresh_device(par);

        soleq{i,j} = equilibrate(par);
        
        Voc_max = 1.2;
        num_points = 281; 
        while Voc_max >= 1.05
            try
                solCV{i, j} = doCV(soleq{i, j}.ion, illumination, -0.2, Voc_max, -0.2, 25e-3, 1, num_points);
                error_log(i,j) = 0;
                results{i,j} = CVstats(solCV{i, j});
                Voc_max = 0;
            catch
                if Voc_max > 1.05
                    warning("Ionic JV solution failed, reducing Vmax by 0.03 V")
                    Voc_max = Voc_max - 0.03;
                    num_points = num_points - 6;
                elseif Voc_max == 1.05
                    warning("Ionic JV solution failed.")
                    error_log(i,j) = 1;
                    results{i,j} = 0;
                end
            end
        end
    end
end

toc

%% Plot results 
Stats_array = zeros(n_k_rec, n_HOMOs, 4);
for i = 1:n_k_rec
    for j = 1:n_HOMOs
        try
            Stats_array(i,j,1) = 1e3*results{i,j}.Jsc_f;
            Stats_array(i,j,2) = results{i,j}.Voc_f;
            Stats_array(i,j,3) = results{i,j}.FF_f;
            Stats_array(i,j,4) = results{i,j}.efficiency_f;
        catch
            warning('No Stats')
            Stats_array(i,j,:) = 0;
        end
    end
end 

%% Run JV for C60 device to get the Jsc value 
parC60 = pc('Input_files/SAM_MAFACsPbIBr_C60.csv');
eqm_QJV_C60 = equilibrate(parC60);
CV_sol_C60 = doCV(eqm_QJV_C60.ion, illumination, -0.2, 1.2, -0.2, 25e-3, 1, 281);
stats_C60 = CVstats(CV_sol_C60);

%%
figure('Name', 'Jsc vs Energy Offsets vs k_rec', 'Position', [50 50 1000 1500])
num = 2;
labels = ["J_{SC} (mA cm^{-2})", "V_{OC} (V)", "FF", "PCE (%)"];
LegendLoc = ["northeast", "southwest", "southeast", "southeast"];
lims = [[-25 -15]; [0.77 1.17]; [0.1, 0.85]; [1 21]];
colormap(flipud('parula'))

box on 
contourf(Donor_HOMO+5.5, k_rec, Stats_array(:,:,num), 'LineStyle', 'none')
hold on
contour(Donor_HOMO+5.5, k_rec, Stats_array(:,:,num), [-1e3*stats_C60.Jsc_f -1e3*stats_C60.Jsc_f], 'color', 'black', 'LineWidth', 3)
hold off
c = colorbar;
c.Label.String = labels(num);
c.Label.FontSize = 30;

set(gca, 'Fontsize', 25, 'YScale', 'log')
xlabel('Energetic Barrier to Hole Injection (eV)', 'FontSize', 30)
ylabel('BHJ k_{2} (cm^{3} s^{-1})', 'FontSize', 30)
xticks([0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4])
xticklabels({'0.10', '0.15', '0.20', '0.25', '0.30', '0.35', '0.40'})
xlim([0.1, 0.4])
ylim([1e-11, 1e-10])

%% Save results and solutions
save_file = 0;
if save_file == 1
    filename = 'PeroBHJ_QuasiTandem_DonorHOMO_BHJkrec_ParameterSweep_CsFAMA_HigherRes.mat';
    save(filename, 'solCV')
end
