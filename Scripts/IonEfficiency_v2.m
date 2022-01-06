tic
%% Define parameter space
%Rows are either recombination velocities or tau_SRH, columns are Ebi values
%Use these vaues for tau_SRH
sigma_values = logspace(-9,-4,6);
%Use these vaues for surface recombination 
%Sigma_values= logspace(-1,5,7);
Ebi_values = zeros(1,14);
Ebi_values(2:end) = linspace(0,0.6,13);
Ebi_values(2) = 0.01;
Nsigma = length(sigma_values); 
NEbi = length(Ebi_values);
params = cell(Nsigma, NEbi);

for i=1:Nsigma
    for j=1:NEbi
            params{i,j} = [sigma_values(i), Ebi_values(j)];
    end
end

error_log = zeros(Nsigma, NEbi);
soleq = cell(Nsigma, NEbi);
solCV_ion = cell(Nsigma, NEbi);
solCV_el = cell(Nsigma, NEbi);
ion_results = cell(Nsigma, NEbi);
el_results = cell(Nsigma, NEbi);
%% Do (many) JV sweeps
par=pc('Input_files/1_layer_test.csv');
for i=1:Nsigma
    for j=13:14%NEbi
        disp(["minority carrier SRV = ", num2str(sigma_values(i)), "cm s-1"])
        disp(["Built in potential, Vbi = ", num2str(2*Ebi_values(j)), " V"])
        par.Phi_left = par.EF0 - params{i,j}(2);
        par.Phi_right = par.EF0 + params{i,j}(2);
        par.taun = params{i,j}(1);
        par.taup = params{i,j}(1);
        %par.sn_l = params{i,j}(1);
        %par.sp_r = params{i,j}(1);
        par = refresh_device(par);

        soleq{i,j} = equilibrate(par);
        try
            if params{i,j}(2) <= 0.35
                solCV_ion{i, j} = doCV(soleq{i, j}.ion, 1, -0.1, 1, -0.1, 1e-4, 1, 241);
                solCV_el{i, j} = doCV(soleq{i, j}.el, 1, -0.1, 1, -0.1, 1e-4, 1, 241);
            elseif params{i,j}(2) > 0.35
                solCV_ion{i, j} = doCV(soleq{i, j}.ion, 1, -0.2, 1.3, -0.2, 1e-4, 1, 241);
                solCV_el{i, j} = doCV(soleq{i, j}.el, 1, -0.2, 1.3, -0.2, 1e-4, 1, 241);
            end
            error_log(i,j) = 0;
            ion_results{i,j} = CVstats(solCV_ion{i, j});
            el_results{i,j} = CVstats(solCV_el{i, j});
        catch
            warning("JV solution failed, try reducing Vmax")
            error_log(i,j) = 1;
        end
    end
end

toc
%% Calculate figure of merit
PCE_ratio = zeros(Nsigma, NEbi);
for i=1:Nsigma
    for j=1:NEbi
        try
            PCE_ratio(i,j) = el_results{i,j}.efficiency_f/ion_results{i,j}.efficiency_f;
        catch
            warning('CVstats unsucessful. Assigning a PCE value of 0.');
            PCE_ratio(i,j) = NaN;
        end
    end
end

%% Plot results 
figure(1)
contourf(2*Ebi_values, log10(sigma_values(1:end)), log10(PCE_ratio(1:end,:)), 25, 'LineWidth', 0.1)
xlabel('V_{BI} (V)')
ylabel('log_{10}(\tau_{SRH} / s)')
%ylabel('log_{10}(s_{surf} /cms^{-1})')
c = colorbar;
c.Label.String = 'log_{10}(PCE_{el} / PCE_{ion})';

%% Voc vs V_bi plot
V_ratio = zeros(Nsigma,NEbi-1);
for i = 1:Nsigma
    for j = 2:NEbi
        V_ratio(i,j) = ion_results{i,j}.Voc_f/(2*Ebi_values(j));
    end
end
figure(2)
contourf(2*Ebi_values(1:end), log10(sigma_values), log10(V_ratio), 12, 'LineWidth', 1)
xlabel('V_{BI} (V)')
%ylabel('log_{10}(s_{surf} /cms^{-1})')
ylabel('log_{10}(\tau_{SRH} / s)')
c = colorbar;
c.Label.String = 'log_{10}(V_{OC,ion} / V_{BI})';