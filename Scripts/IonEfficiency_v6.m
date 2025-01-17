%Use this file to symmetrically sweep HOMO/LUMO offsets vs ion
%concentration for a device with doped or undoped interlayers

%There have been some choices made here
%Doping is fixed s.t. Fermi Level offset is alway 0.1 eV from the relevant
%band edge
%Work function of the electrodes is handled as follows
%Either the Vbi is 1.1 eV, with an equal offset from the ETL and HTL CB/VB,
%respectively
%Or, if the offset is so large that a Vbi of 1.1 eV would mean the electode
%work function lying within 0.1 eV of the CB/VB edge, the work functions
%are redefined so that they lie 0.1 eV from the CB/VB edge 

%Used this file to produce v4 of the simulation results 

%TURN SAVE OFF TO START OFF WITH (final cell)

tic
%% Define parameter space
%Choose to use doped or undoped TLs
doped = 1;
high_performance = 0;
n_values = 7;
Delta_TL = linspace(0, 0.3, n_values);
% n_values = 2;
% Delta_TL = [0 0.3];
Symmetric_offset = 0;
%Fix the offset for the ETL or HTL
Fix_ETL = 1;
%Energetic offset between the perovskite and TL for the TL with fixed
%energetics 
Fixed_offset = 0.15;
%This is a bit of a hack, but if the offfset is exactly 0, the surface
%recombination error becomes huge for reasons I do not fully understand...
%The saga continues - only seems to matter for the HTL, not the ETL...
Delta_TL(1) = 1e-3;
Ion_Conc = [1e15 5e15 1e16 5e16 1e17 5e17 1e18 0];
% Ion_Conc = [1e18 0];
n_ion_concs = length(Ion_Conc);

%Rows are the Ion Concentrations
%Columns are the TL Energetic Offsets
params = cell(n_ion_concs, n_values);

for i=1:n_ion_concs
    for j=1:n_values
            params{i,j} = [Ion_Conc(i), Delta_TL(j)];
    end
end

%% Set up structures for storing the results
error_log = zeros(n_ion_concs, n_values);
soleq = cell(n_ion_concs, n_values);
solCV = cell(n_ion_concs, n_values);
results = cell(n_ion_concs, n_values);

%% Do (many) JV sweeps
%Select the correct input file for doped or undoped cases 
if doped == 1
    if high_performance == 0
        par=pc('Input_files/EnergyOffsetSweepParameters_v5_doped.csv');
    elseif high_performance == 1
        par=pc('Input_files/EnergyOffsetSweepParameters_v5_doped_higherPCE.csv');
    end
elseif doped == 0
    if high_performance == 0
        par=pc('Input_files/EnergyOffsetSweepParameters_v5_undoped.csv');
    elseif high_performance == 1
        par=pc('Input_files/EnergyOffsetSweepParameters_v5_undoped_higherPCE.csv');
    end
end 

%Set the illumination for the JV sweeps and the max voltage to go up to in
%the JV sweeps
if high_performance == 0
    illumination = 1;
    max_val = 1.2;
    phi_L = -5.15;
    phi_R = -4.05;
elseif high_performance == 1
    illumination = 1.2;
    max_val = 1.3;
    phi_L = -5.2;
    phi_R = -4.1;
end

%Piers version - set this to 1 to ensure that the WFs of the electrodes is
%always equal to the Fermi level of the contacts - stops there being a
%Schottky barrier at the TL/metal interface, which negatively affects
%device performance. 
%Less similar to the experimental situation maybe, as there the same
%electrode materials are used for both interlayers
Piers_version = 0;
%Only makes sense to do this for the doped case so put this is as a
%protection in case you forget to change it
if doped == 0
    Piers_version = 0;
end

%Reset the electrode work functions in each loop to be safe as they are
%changed for the cases where E_LUMO (E_HOMO) is far below (above) the CB
%(VB)
for i = 1:n_ion_concs
    for j = 1:n_values
        disp(["DeltaE_HOMO = ", num2str(Delta_TL(j)), " eV"])
        if i < n_ion_concs 
            disp(["Ion_Conc = ", num2str(Ion_Conc(i)), " cm{^-3}"])
        elseif i == n_ion_concs 
            disp("No Mobile Ions")
        end

        %HTL Energetics
        if Symmetric_offset == 1 || (Symmetric_offset == 0 && Fix_ETL == 1)
            par.Phi_left = phi_L;
            par.Phi_IP(1) = par.Phi_IP(3) + params{i,j}(2);
            par.Phi_EA(1) = par.Phi_IP(1) + 2.5;
            par.Et(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
            if doped == 0
                par.EF0(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
            elseif doped == 1
                par.EF0(1) = par.Phi_IP(1) + 0.1;
            end
            if Piers_version == 0
                if par.Phi_left < par.Phi_IP(1) + 0.1
                    par.Phi_left = par.Phi_IP(1) + 0.1;
                end
            elseif Piers_version == 1
                par.Phi_left = par.Phi_IP(1) + 0.1;
            end
        elseif Symmetric_offset == 0 && Fix_ETL == 0
            par.Phi_left = phi_L;
            par.Phi_IP(1) = par.Phi_IP(3) + Fixed_offset;
            par.Phi_EA(1) = par.Phi_IP(1) + 2.5;
            par.Et(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
            if doped == 0
                par.EF0(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
            elseif doped == 1
                par.EF0(1) = par.Phi_IP(1) + 0.1;
            end
            if Piers_version == 0
                if par.Phi_left < par.Phi_IP(1) + 0.1
                    par.Phi_left = par.Phi_IP(1) + 0.1;
                end
            elseif Piers_version == 1
                par.Phi_left = par.Phi_IP(1) + 0.1;
            end
        end

        %ETL Energetics
        %Need to use opposite sign at ETL to keep energy offsets symmetric
        if Symmetric_offset == 1 || (Symmetric_offset == 0 && Fix_ETL == 0)
            par.Phi_right = phi_R;
            par.Phi_EA(5) = par.Phi_EA(3) - params{i,j}(2);
            par.Phi_IP(5) = par.Phi_EA(5) - 2.5;
            par.Et(5) = (par.Phi_IP(5) + par.Phi_EA(5))/2;
            if doped == 0
                par.EF0(5) = (par.Phi_IP(5) + par.Phi_EA(5))/2;
            elseif doped == 1
                par.EF0(5) = par.Phi_EA(5) - 0.1;
            end
            if Piers_version == 0
                if par.Phi_right > par.Phi_EA(5) - 0.1
                    par.Phi_right = par.Phi_EA(5) - 0.1;
                end
            elseif Piers_version == 1
                par.Phi_right = par.Phi_EA(5) - 0.1;
            end
        elseif Symmetric_offset == 0 && Fix_ETL == 1
            par.Phi_right = phi_R;
            par.Phi_EA(5) = par.Phi_EA(3) - Fixed_offset;
            par.Phi_IP(5) = par.Phi_EA(5) - 2.5;
            par.Et(5) = (par.Phi_IP(5) + par.Phi_EA(5))/2;
            if doped == 0
                par.EF0(5) = (par.Phi_IP(5) + par.Phi_EA(5))/2;
            elseif doped == 1
                par.EF0(5) = par.Phi_EA(5) - 0.1;
            end
            if Piers_version == 0
                if par.Phi_right > par.Phi_EA(5) - 0.1
                    par.Phi_right = par.Phi_EA(5) - 0.1;
                end
            elseif Piers_version == 1
                par.Phi_right = par.Phi_EA(5) - 0.1;
            end
        end

        %ion conc
        if i ~= n_ion_concs
            par.Ncat(:) = params{i,j}(1);
            par.Nani(:) = params{i,j}(1);
        end 
        
        %Do this as it seesm to reduce the discrepency between surface
        %volumetric surace recombination model and the abrupt interface one
        %But also made solution less stable? Maybe better to tinker with
        %this on the one which varies surface recombination...
        %par.frac_vsr_zone = 0.05;
        par = refresh_device(par);

        soleq{i,j} = equilibrate(par);
        
        %electron only scan
        if i == n_ion_concs 
            Fermi_offset = par.EF0(5) - par.EF0(1);
            if Fermi_offset > max_val
                Voc_max = Fermi_offset + 0.05;
            else
                Voc_max = max_val;
            end
            num_points = (2*100*(Voc_max+0.2))+1;
            while Voc_max >= 1.05
                try            
                    solCV{i, j} = doCV(soleq{i, j}.el, illumination, -0.2, Voc_max, -0.2, 1e-4, 1, num_points);           
                    error_log(i,j) = 0;
                    results{i,j} = CVstats(solCV{i, j});
                    Voc_max = 0;                
                catch
                    if Voc_max > 1.05
                        warning("Electronic-only JV solution failed, reducing Vmax by 0.03 V")
                        Voc_max = Voc_max - 0.03;
                        num_points = num_points - 6;
                    elseif Voc_max == 1.05
                        warning("Electronic-only JV solution failed.")
                        error_log(i,j) = 1;
                        results{i,j} = 0;
                    end
                end
            end
        
        else
            Fermi_offset = par.EF0(5) - par.EF0(1);
            if Fermi_offset > max_val
                Voc_max = Fermi_offset + 0.05;
            else
                Voc_max = max_val;
            end
            num_points = (2*100*(Voc_max+0.2))+1; 
            while Voc_max >= 1.05
                try
                    solCV{i, j} = doCV(soleq{i, j}.ion, illumination, -0.2, Voc_max, -0.2, 1e-4, 1, num_points);
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
end

toc

%% Plot results 
plot_JVs = 1;
Stats_array = zeros(n_ion_concs, n_values, 5);
if plot_JVs == 1
    J_result = cell(2, n_values);
%     J_srh_result = cell(2, n_values);
%     J_vsr_result = cell(2, n_values);
    V_plot = cell(2, n_values);
end
e = solCV{1,1}.par.e;

for i = 1:n_ion_concs
    for j = 1:n_values
        try
            Stats_array(i,j,1) = 1e3*results{i,j}.Jsc_f;
            Stats_array(i,j,2) = results{i,j}.Voc_f;
            Stats_array(i,j,3) = results{i,j}.FF_f;
            Stats_array(i,j,4) = results{i,j}.efficiency_f;
            J_result{i,j} = dfana.calcJ(solCV{i,j}).tot(:,1);
            V_plot{i,j} = dfana.calcVapp(solCV{i,j});
%             x = solCV{i,j}.par.x_sub;
%             loss_currents = dfana.calcr(solCV{i,j},'sub');
%             Vapp = dfana.calcVapp(solCV{i,j});
%             end_value = cast((length(Vapp)-1)/2, 'int32');
%             J_srh = e*trapz(x, loss_currents.srh, 2)';
%             J_vsr = e*trapz(x, loss_currents.vsr, 2)';
%             J_srh_OC = interp1(Vapp(1:end_value), J_srh(1:end_value), Stats_array(i,j,2));
%             J_vsr_OC = interp1(Vapp(1:end_value), J_vsr(1:end_value), Stats_array(i,j,2));
%             if J_srh_OC > J_vsr_OC
%                 Stats_array(i,j,5) = 0;
%             elseif J_srh_OC < J_vsr_OC
%                 Stats_array(i,j,5) = 1;
%             end
%             if plot_JVs == 1
%                 if i == 7
%                     V_plot{1,j} = Vapp;
%                     J_srh_result{1,j} = J_srh;
%                     J_vsr_result{1,j} = J_vsr;
%                 elseif i == 8
%                     V_plot{2,j} = Vapp;
%                     J_srh_result{2,j} = J_srh;
%                     J_vsr_result{2,j} = J_vsr;
%                 end
%             end
        catch
            warning('No Stats')
            Stats_array(i,j,:) = 0;
        end
    end
end 

%%
figure('Name', 'JV Parameter vs Energy Offsets vs Ion Conc', 'Position', [50 50 800 800])
Colours = parula(n_ion_concs-1);
num = 2;
labels = ["J_{SC} (mA cm^{-2})", "V_{OC} (V)", "FF", "PCE (%)"];
LegendLoc = ["northeast", "southwest", "southeast", "northeast"];
if doped == 0
%     lims = [[-23 -15]; [0.77 1.24]; [0.5, 0.9]; [10 23]];
    lims = [[-23 -15]; [0.77 1.24]; [0.5, 0.9]; [16 27]];
elseif doped == 1
%     lims = [[-24 -15]; [0.77 1.24]; [0.5, 0.9]; [10 23]];
    lims = [[-24 -15]; [0.77 1.24]; [0.5, 0.9]; [18 27]];
end
box on 
for i = 1:n_ion_concs
    hold on
    if i == n_ion_concs
        plot(Delta_TL, Stats_array(n_ion_concs,:,5).*Stats_array(n_ion_concs,:,num), 'marker', 'x', 'Color', 'black', 'LineStyle', 'none', 'MarkerSize', 10, 'HandleVisibility', 'Off')
        plot(Delta_TL, (1-Stats_array(n_ion_concs,:,5)).*Stats_array(n_ion_concs,:,num), 'marker', 'o', 'Color', 'black', 'LineStyle', 'none', 'MarkerSize', 10, 'HandleVisibility', 'Off')
        plot(Delta_TL, Stats_array(n_ion_concs,:,num), 'marker', 'none', 'Color', 'black')
    else
%         plot(Delta_TL, Stats_array(i-1,:,5).*Stats_array(i-1,:,num), 'marker', 'x', 'Color', Colours(i-1,:), 'LineStyle', 'none', 'MarkerSize', 10, 'HandleVisibility', 'Off')
%         plot(Delta_TL, (1-Stats_array(i-1,:,5)).*Stats_array(i-1,:,num), 'marker', 'o', 'Color', Colours(i-1,:), 'LineStyle', 'none', 'MarkerSize', 10, 'HandleVisibility', 'Off')
        plot(Delta_TL, Stats_array(i,:,num), 'marker', 'none', 'Color', Colours(i,:))
    end
end
set(gca, 'Fontsize', 25)
xlabel('Transport Layer Energetic Offset (eV)', 'FontSize', 30)
ylabel(labels(num), 'FontSize', 30)
xlim([0, 0.3])
xticks([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3])
xticklabels({'0.00', '0.05', '0.10', '0.15', '0.20', '0.25', '0.30'})
ylim(lims(num,:))
%legend({'1e15', '5e15', '1e16', '5e16', '1e17', '5e17', '1e18', 'No Ions'}, 'Location', LegendLoc(num), 'FontSize', 25, 'NumColumns', 2)
title(legend, 'Ion Concentration (cm^{-3})', 'FontSize', 25)

%% Plot JV curves w/w.o. ions
plot_JVs = 1;
if plot_JVs == 1

    figure('Name', 'JVPlot', 'Position', [100 100 800 800])
    %Colours = flip(parula(n_values));
    Colours = {[0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0 0.4470 0.7410], [0.4940 0.1840 0.5560]};
    %Ion_Conc = [1e15 5e15 1e16 5e16 1e17 5e17 1e18 0];
    %Set which ion concentration to plot for
    %Have coppied the array above so you can see which nuber is the right one
    %easily
    num_Ion_Conc = 2;
    
    count = 1;
    for j = [1 3 5 7]
        v = dfana.calcVapp(solCV{num_Ion_Conc, j});
        J = dfana.calcJ(solCV{num_Ion_Conc, j}).tot(:,1);
        
        hold on
        xline(0, 'black', 'HandleVisibility', 'off')
        yline(0, 'black', 'HandleVisibility', 'off')
        plot(v(:), J(:)*1000, 'color', Colours{count}, 'LineWidth', 3) 
        count = count + 1;
        hold off
    
    end
    
    box on 
    set(gca, 'FontSize', 25)
    xlim([-0.15, 1.2])
    ylim([-25,5])

    legend({'  0.0', '  0.1', '  0.2', '  0.3'}, 'Location', 'northwest', 'FontSize', 25, 'NumColumns', 2)
    title(legend, ['Transport Layer' newline 'Energetic Offset (eV)'], 'FontSize', 25)
    xlabel('Voltage(V)', 'FontSize', 30)
    ylabel('Current Density (mAcm^{-2})', 'FontSize', 30)
    ax1 = gcf;
    
end

%% Plot J_srh and _sr for high ions vs no ions as a function of offset
plot_JVs = 1;
if plot_JVs == 1
    figure('Name', 'JVPlot', 'Position', [100 100 800 800])
    Colours = parula(n_values);
    
    for j = 1:n_values
        
        subplot(1,2,1)
        
         xline(0, 'black', 'HandleVisibility', 'off')
         yline(0, 'black', 'HandleVisibility', 'off')
        if j == 1
            plot(V_plot{1,j}, J_vsr_result{1,j}*1000, 'color', Colours(j,:), 'LineWidth', 3) 
            hold on
            plot(V_plot{2,j}, J_vsr_result{2,j}*1000, 'color', Colours(j,:), 'LineWidth', 3, 'LineStyle', '--', 'HandleVisibility', 'off') 
            hold off
        else
            hold on
            plot(V_plot{1,j}, J_vsr_result{1,j}*1000, 'color', Colours(j,:), 'LineWidth', 3) 
            plot(V_plot{2,j}, J_vsr_result{2,j}*1000, 'color', Colours(j,:), 'LineWidth', 3, 'LineStyle', '--', 'HandleVisibility', 'off') 
            hold off
        end
            
        box on
    
        set(gca, 'FontSize', 25)
        xlim([-0.15, 1.2])
        ylim([-5, 25])
        
        legend({'0.00', '0.05', '0.10', '0.15', '0.20', '0.25', '0.30'}, 'Location', 'southwest', 'FontSize', 25)
        title(legend, 'Transport Layer Offset (eV)', 'FontSize', 25)
        
        xlabel('Voltage(V)', 'FontSize', 30)
        ylabel('Current Density (mAcm^{-2})', 'FontSize', 30)
    
    end
    
    for j = 1:n_values
        
        subplot(1,2,2)
    
         xline(0, 'black', 'HandleVisibility', 'off')
         yline(0, 'black', 'HandleVisibility', 'off')
    
        if j == 1
            plot(V_plot{1,j}, J_srh_result{1,j}*1000, 'color', Colours(j,:), 'LineWidth', 3) 
            hold on
            plot(V_plot{2,j}, J_srh_result{2,j}*1000, 'color', Colours(j,:), 'LineWidth', 3, 'LineStyle', '--', 'HandleVisibility', 'off') 
            hold off
        else
            hold on
            plot(V_plot{1,j}, J_srh_result{1,j}*1000, 'color', Colours(j,:), 'LineWidth', 3) 
            plot(V_plot{2,j}, J_srh_result{2,j}*1000, 'color', Colours(j,:), 'LineWidth', 3, 'LineStyle', '--', 'HandleVisibility', 'off') 
            hold off
        end
        box on 
    
        set(gca, 'FontSize', 25)
        xlim([-0.15, 1.2])
        ylim([-5, 25])
        
        legend({'0.00', '0.05', '0.10', '0.15', '0.20', '0.25', '0.30'}, 'Location', 'southwest', 'FontSize', 25)
        title(legend, 'Transport Layer Offset (eV)', 'FontSize', 25)
        
        xlabel('Voltage(V)', 'FontSize', 30)
        ylabel('Current Density (mAcm^{-2})', 'FontSize', 30)
    
    end
end

%% Save results and solutions
save_file = 1;
if save_file == 1
    if doped == 0
        filename = 'DeltaE_v5_undoped_SAMComparison.mat';
    elseif doped == 1
        filename = 'DeltaE_v5_doped_HigherPerforamnce.mat';
    end 
    save(filename, 'results', 'solCV')
end
