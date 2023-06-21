%Use this file to symmetrically sweep HOMO/LUMO offsets vs surface
%recombination velocity

%TURN SAVE OFF TO START OFF WITH (final cell)

tic
%% Define parameter space
%Rows are the surface recombination velocities
%Columns are the TL Energetic Offsets

n_values = 10;
Delta_TL = linspace(-0.15, 0.3, n_values);
v_sr = [0.5, 5, 50, 500, 5000];
n_v_sr  = length(v_sr);

params = cell(n_v_sr, n_values);

for i=1:n_v_sr
    for j=1:n_values
            params{i,j} = [v_sr(i), Delta_TL(j)];
    end
end

%%
soleq = cell(n_v_sr, n_values);

error_log = zeros(n_v_sr, n_values);
solCV = cell(n_v_sr, n_values);
results = cell(n_v_sr, n_values);

error_log_el = zeros(n_v_sr, n_values);
solCV_el = cell(n_v_sr, n_values);
results_el = cell(n_v_sr, n_values);

%% Do (many) JV sweeps
%Remeber to update the work function values if you change these parameters
%between files 
par=pc('Input_files/EnergyOffsetSweepParameters_v3.csv');
illumination = 1.1;

%Reset the electrode work functions in each loop to be safe as they are
%changed for the cases where E_LUMO (E_HOMO) is far below (above) the CB
%(VB)
for i = 1:n_v_sr
    for j = 1:n_values
        disp(["DeltaE_HOMO = ", num2str(Delta_TL(j)), " eV"])
        disp(["v_sr = ", num2str(v_sr(i)), " cm s{^-1}"])
        %HTL Energetics
        par.Phi_left = -5.15;
        par.Phi_IP(1) = par.Phi_IP(3) + params{i,j}(2);
        par.Phi_EA(1) = par.Phi_IP(1) + 2.5;
        par.EF0(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
        par.Et(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
        if par.Phi_left < par.Phi_IP(1)
            par.Phi_left = par.Phi_IP(1);
        end
        %ETL Energetics
        %Need to use opposite sign at ETL to keep energy offsets symmetric
        par.Phi_right = -4.05;
        par.Phi_EA(5) = par.Phi_EA(3) - params{i,j}(2);
        par.Phi_IP(5) = par.Phi_EA(5) - 2.5;
        par.EF0(5) = (par.Phi_IP(5)+par.Phi_EA(5))/2;
        par.Et(5) = (par.Phi_IP(5)+par.Phi_EA(5))/2;
        if par.Phi_right > par.Phi_EA(5)
            par.Phi_right = par.Phi_EA(5);
        end
        %v_sr
        par.sn(2) = params{i,j}(1);
        par.sn(4) = params{i,j}(1);
        par.sp(2) = params{i,j}(1);
        par.sp(2) = params{i,j}(1);

        par = refresh_device(par);

        soleq{i,j} = equilibrate(par);
        
        %electron only scan
        Voc_max = 1.2;
        num_points = 281;
        while Voc_max >= 1.05
            try            
                solCV_el{i, j} = doCV(soleq{i, j}.el, illumination, -0.2, Voc_max, -0.2, 1e-4, 1, num_points);           
                error_log_el(i,j) = 0;
                results_el{i,j} = CVstats(solCV_el{i, j});
                Voc_max = 0;                
            catch
                if Voc_max > 1.05
                warning("Electronic-only JV solution failed, reducing Vmax by 0.03 V")
                Voc_max = Voc_max - 0.03;
                num_points = num_points - 6;
                elseif Voc_max == 1.05
                    warning("Electronic-only JV solution failed.")
                    error_log_el(i,j) = 1;
                    results_el{i,j} = 0;
                end
            end
        end
        
        Voc_max = 1.2;
        num_points = 281; 
        while Voc_max >= 1.05
            try
                solCV{i, j} = doCV(soleq{i, j}.ion, illumination, -0.2, Voc_max, -0.2, 1e-4, 1, num_points);
                error_log(i,j) = 0;
                results{i,j} = CVstats(solCV{i, j});
                Voc_max = 0;
            catch
                if Voc_max > 1.05
                    warning("Ionic JV solution failed, reducing Vmax by 0.05 V")
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
Stats_array = zeros(n_v_sr, n_values, 4);
for i = 1:n_v_sr
    for j = 1:n_values
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

Stats_array_el = zeros(n_v_sr, n_values, 4);
for i = 1:n_v_sr
    for j = 1:n_values
        try
            Stats_array_el(i,j,1) = 1e3*results_el{i,j}.Jsc_f;
            Stats_array_el(i,j,2) = results_el{i,j}.Voc_f;
            Stats_array_el(i,j,3) = results_el{i,j}.FF_f;
            Stats_array_el(i,j,4) = results_el{i,j}.efficiency_f;
        catch
            warning('No Stats')
            Stats_array(i,j,:) = 0;
        end
    end
end 
%%
figure('Name', 'PCE vs Energy Offsets vs Ion Conc', 'Position', [50 50 1000 1200])
Colours = parula(n_v_sr);
num = 2;
labels = ["J_{SC} (mA cm^{-2})", "V_{OC} (V)", "FF", "PCE (%)"];
LegendLoc = ["northeast", "southwest", "southeast", "southeast"];
lims = [[-24 -15]; [0.77 1.24]; [0.1, 0.85]; [1 21]];
box on 
for i = 1:n_v_sr
    hold on
    plot(Delta_TL, Stats_array(i,:,num), 'marker', 'x', 'Color', Colours(i,:))
    plot(Delta_TL, Stats_array_el(i,:,num), 'marker', 'o', 'Color', Colours(i,:), 'LineStyle', ':')
end
set(gca, 'Fontsize', 25)
xlabel('Transport Layer Energetic Offset (eV)', 'FontSize', 30)
ylabel(labels(num), 'FontSize', 30)
xlim([-0.15, 0.3])
ylim(lims(num,:))
legend({'0.5', '', '5', '', '50', '', '500', '', '5000', ''}, 'Location', LegendLoc(num), 'FontSize', 25, 'NumColumns', 2)
title(legend, 'Surface Recombination Velocity (cm s^{-1})', 'FontSize', 25)

%% Save results and solutions
save_file = 0;
if save_file == 1
    filename = 'tld_symmetric_DeltaEHOMO_vs_DeltaELUMO_v2.mat';
    save(filename, 'el_results', 'results', 'solCV_el', 'solCV')
end
