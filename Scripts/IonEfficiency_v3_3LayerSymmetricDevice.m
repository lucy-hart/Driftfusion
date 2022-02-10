tic
%% Define parameter space
%Rows are variable parameter, columns are the ion concentrations

%Use these vaues for surface recombination 
sigma_values= logspace(-1,5,7);
IonConc_values = [1e15, 1e16, 1e17, 1e19];
Nsigma = length(sigma_values); 
NIonConc = length(IonConc_values);
params = cell(Nsigma, NIonConc);

for i=1:Nsigma
    for j=1:NIonConc
            params{i,j} = [sigma_values(i), IonConc_values(j)];
    end
end

error_log = zeros(Nsigma, NIonConc);
soleq = cell(Nsigma, NIonConc);
solCV_ion = cell(Nsigma, NIonConc);
solCV_el = cell(Nsigma, NIonConc);
ion_results = cell(Nsigma, NIonConc);
el_results = cell(Nsigma, NIonConc);

%% Do (many) JV sweeps
par=pc('Input_files/3_layer_test_symmetric.csv');
for i=1:Nsigma
    for j=1:NIonConc
        disp(["minority carrier SRV = ", num2str(sigma_values(i)), "cm s-1"])
        disp(["Ion conc, N_cat = ", num2str(IonConc_values(j)), " cm-3"])
        par.Ncat(:) = params{i,j}(2);
        par.Nani(:) = params{i,j}(2);
        par.sn(2) = params{i,j}(1);
        par.sp(4) = params{i,j}(1);
        par = refresh_device(par);

        soleq{i,j} = equilibrate(par);
        try
            solCV_ion{i, j} = doCV(soleq{i, j}.ion, 1, -0.1, 1.4, -0.1, 1e-4, 1, 241);
            solCV_el{i, j} = doCV(soleq{i, j}.el, 1, -0.1, 1.4, -0.1, 1e-4, 1, 241);
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
PCE_ratio = zeros(Nsigma, NIonConc);
for i=1:Nsigma
    for j=1:NIonConc
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
contourf(log10(IonConc_values), log10(sigma_values(1:end)), log10(PCE_ratio(1:end,:)), 25, 'LineWidth', 0.1)
xlabel('log_{10}(Mobile Cation Concentration /cm^{-3})')
ylabel('log_{10}(s_{surf} /cms^{-1})')
c = colorbar;
c.Label.String = 'log_{10}(PCE_{el} / PCE_{ion})';

%% Voc vs V_bi plot
V_ratio = zeros(Nsigma,NIonConc-1);
for i = 1:Nsigma
    for j = 2:NIonConc
        V_ratio(i,j) = ion_results{i,j}.Voc_f/(2*Ebi_values(j));
    end
end
figure(2)
contourf(2*Ebi_values(1:end), log10(sigma_values), log10(V_ratio), 12, 'LineWidth', 1)
xlabel('log_{10}(Mobile Cation Concentration /cm^{-3})')
ylabel('log_{10}(s_{surf} /cms^{-1})')
c = colorbar;
c.Label.String = 'log_{10}(V_{OC,ion} / V_{BI})';