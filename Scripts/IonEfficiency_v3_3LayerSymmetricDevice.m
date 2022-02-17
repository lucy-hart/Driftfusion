tic
%% Define parameter space
%Rows are variable parameter, columns are the ion concentrations

IonConc_values = [1e15, 1e16, 1e17, 1e18, 1e19];
%Use these vaues for surface recombination 
%y_values= logspace(-1,5,7);
%Use these vaues for tau_SRH
y_values = logspace(-9,-4,6);

N_y = length(y_values); 
NIonConc = length(IonConc_values);
params = cell(N_y, NIonConc);

for i=1:N_y
    for j=1:NIonConc
            params{i,j} = [y_values(i), IonConc_values(j)];
    end
end

error_log = zeros(N_y, NIonConc);
soleq = cell(N_y, NIonConc);
solCV_ion = cell(N_y, NIonConc);
solCV_el = cell(N_y, NIonConc);
ion_results = cell(N_y, NIonConc);
el_results = cell(N_y, NIonConc);

%% Do (many) JV sweeps
par=pc('Input_files/3_layer_test_symmetric.csv');
for i=1:N_y
    for j=1:NIonConc
        disp(["tau_SRH = ", num2str(y_values(i)), " s"])
        %disp(["minority carrier SRV = ", num2str(y_values(i)), "cm s-1"])
        disp(["Ion conc, N_cat = ", num2str(IonConc_values(j)), " cm-3"])
        par.Ncat(:) = params{i,j}(2);
        par.Nani(:) = params{i,j}(2);
        %Remember to change these when you change y variable
        par.taun(3) = params{i,j}(1);
        par.taup(3) = params{i,j}(1);
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
PCE_ratio = zeros(N_y, NIonConc);
for i=1:N_y
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
contourf(log10(IonConc_values), log10(y_values(1:end)), log10(PCE_ratio(1:end,:)), 25, 'LineWidth', 0.1)
xlabel('log_{10}(Mobile Cation Concentration /cm^{-3})')
%ylabel('log_{10}(s_{surf} /cms^{-1})')
ylabel('log_{10}(\tau_{SRH} /s)')
c = colorbar;
c.Label.String = 'log_{10}(PCE_{el} / PCE_{ion})';

%% Save results and solutions

filename = 'tld_symmetric_tSRH_vs_Ncat.mat';
save(filename, 'el_results', 'ion_results', 'solCV_el', 'solCV_ion')
