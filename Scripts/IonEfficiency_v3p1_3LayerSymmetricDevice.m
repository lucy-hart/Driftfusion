%Use this file for energy related things (TL Fermi levels, energy offsets
%and Vbi)

tic
%% Define parameter space
%Rows are variable parameter, columns are the ion concentrations

IonConc_values = [1e15, 1e16, 1e17, 1e18, 1e19];
%Use these vaues for Vbi 
y_values= zeros(1,14);
y_values(2:end) = linspace(0,0.6,13);
y_values(2) = 0.01;

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
%For SRH reset tau values in transport layers too
for i=1:N_y
    for j=1:NIonConc
        disp(["Ion Concentration = ", num2str(IonConc_values(i)), " cm-3"])
        disp(["Built in potential, Vbi = ", num2str(2*y_values(j)), " V"])
        par.Ncat(:) = params{i,j}(2);
        par.Nani(:) = params{i,j}(2);
        %Remember to change these when you change y variable
        par.Phi_left = par.EF0(3) - params{i,j}(1);
        par.Phi_right = par.EF0(3) + params{i,j}(1);
        par = refresh_device(par);

        soleq{i,j} = equilibrate(par);
        try
            if params{i,j}(1) <= 0.35
                solCV_ion{i, j} = doCV(soleq{i, j}.ion, 1, -0.1, 1, -0.1, 1e-4, 1, 241);
                solCV_el{i, j} = doCV(soleq{i, j}.el, 1, -0.1, 1, -0.1, 1e-4, 1, 241);
            elseif params{i,j}(1) > 0.35
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
contourf(log10(IonConc_values), 2*y_values, log10(PCE_ratio(1:end,:)), 25, 'LineWidth', 0.1)
xlabel('log_{10}(Mobile Cation Concentration /cm^{-3})')
%ylabel('log_{10}(s_{surf} /cms^{-1})')
ylabel('log_{10}(V_{BI} /V)')
ylim([-7,-4])
c = colorbar;
c.Label.String = 'log_{10}(PCE_{el} / PCE_{ion})';

%% Save results and solutions

filename = 'tld_symmetric_Vbi_vs_Ncat_DopedTLs.mat';
save(filename, 'el_results', 'ion_results', 'solCV_el', 'solCV_ion')
