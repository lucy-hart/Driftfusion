%Use this file to sweep ion concentration vs surface recombination veocity
%or SRH lieftime 

%Can also choose at start what energetic offset to look at 
%Tend to use a relatively large offset as this seems to be the case in the
%systems which Fraser is looking at 

%Despite being an earlier version than v6, this file is still for use with
%v4 of the parameters i.e., it handles both the doped and undoped cases 

%TURN SAVE OFF TO START OFF WITH (final cell)

tic
%% Define parameter space
%Choose to use doped or undoped TLs and which of v_sr or tau_SRH to vary
doped = 1;
surface = 1;

%Set up the parameters for the ion concentrations
Ion_Conc = [5e15 1e16 5e16 1e17 5e17 1e18 0];
%Ion_Conc = [1e18 0];
n_ion_concs = length(Ion_Conc);

if surface == 1
    v_sr = [1];%logspace(-1,4,11);
    n_recom  = length(v_sr);
    %Rows are the ion concentrations    
    %Columns are the surface recombination velocities
    params = cell(n_ion_concs, n_recom);
    for i=1:n_ion_concs
        for j=1:n_recom
            params{i,j} = [Ion_Conc(i), v_sr(j)];
        end
    end
elseif surface == 0
    tau_SRH = [1e-8 3e-8 5e-8 1e-7 3e-7 5e-7 1e-6];
    n_recom  = length(tau_SRH);
    %Rows are the ion concentrations    
    %Columns are the SRH lifetimes 
    params = cell(n_ion_concs, n_recom);
    for i=1:n_ion_concs
        for j=1:n_recom
            params{i,j} = [Ion_Conc(i), tau_SRH(j)];
        end
    end
end

%%
%Select the correct input file for doped or undoped cases
if doped == 0
    par=pc('Input_files/EnergyOffsetSweepParameters_v5_undoped.csv');
    Voc_max_lim = 1.05;
elseif doped == 1
    par=pc('Input_files/EnergyOffsetSweepParameters_v5_doped.csv');
    Voc_max_lim = 1.05;
elseif doped == 0.5
    par=pc('Input_files/TiO2_MAPI_Spiro_TestSaP.csv');
    Voc_max_lim = 0.80;
end

%% Choose the nergetics of the TLs 
%Default values are FILL THIS IN
%Will use these values if Fiddle_with_Energetics is 0
Fiddle_with_Energetics = 1;

if Fiddle_with_Energetics == 1
    %Choose the offsets for the system
    %Positive offset for DHOMO means TL VB lies above the perovskite VB
    %Negative offset for DLUMO means TL CB lies below the perovskite CB
    DHOMO = 0.2;
    DLUMO = -0.2;

    %HTL Energetics
    par.Phi_left = -5.15;
    par.Phi_IP(1) = par.Phi_IP(3) + DHOMO;
    par.Phi_EA(1) = par.Phi_IP(1) + 2.5;
    par.Et(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
    if doped == 0
        par.EF0(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
    elseif doped == 1
        par.EF0(1) = par.Phi_IP(1) + 0.1;
    end 
    if par.Phi_left < par.Phi_IP(1) + 0.01
        par.Phi_left = par.Phi_IP(1) + 0.01;
    end

    %ETL Energetics
    par.Phi_right = -4.05;
    par.Phi_EA(5) = par.Phi_EA(3) + DLUMO;
    par.Phi_IP(5) = par.Phi_EA(5) - 2.5;
    par.Et(5) = (par.Phi_IP(5)+par.Phi_EA(5))/2;
    if doped == 0
         par.EF0(5) = (par.Phi_IP(5)+par.Phi_EA(5))/2;
    elseif doped == 1
        par.EF0(5) = par.Phi_EA(5) - 0.1;
    end
    if par.Phi_right > par.Phi_EA(5) - 0.01
        par.Phi_right = par.Phi_EA(5) - 0.01;
    end

    par = refresh_device(par);

end

%% Set up structures for storing the results
error_log = zeros(n_ion_concs, n_recom);
soleq = cell(n_ion_concs, n_recom);
solCV = cell(n_ion_concs, n_recom);
results = cell(n_ion_concs, n_recom);

%% Do (many) JV sweeps

%Set the illumnation for the JV sweeps 
illumination = 1;

%Reset the electrode work functions in each loop to be safe as they are
%changed for the cases where E_LUMO (E_HOMO) is far below (above) the CB
%(VB)
for i = 1:n_ion_concs
    for j = 1:n_recom
        if surface == 0
            disp(["tau_SRH = ", num2str(tau_SRH(j)), " s"])
        elseif surface == 1
            disp(["v_sr = ", num2str(v_sr(j)), " cm{^-1}"])
        end 
        if i < n_ion_concs 
            disp(["Ion_Conc = ", num2str(Ion_Conc(i)), " cm{^-3}"])
        elseif i == n_ion_concs 
            disp("No Mobile Ions")
        end

        if surface == 0
            par.taun(3) = params{i,j}(2);
            par.taup(3) = params{i,j}(2);
        %NB: set electron and hole surface recombination velocities to be the
        %same as I realised the 'majority' carrier surface recombination
        %velocity can actualy affect the results quite a lot as the lack of
        %field means that there isn't really a well defined
        %maority/mminority carrier at the interface
        %Choose to set them the same since the trap states I use are midgap
        %and so to reason to assume they preferentially trap one carrier
        %type or the other
        %This is an assumption
        %It may very well be wrong
        %Need to remeber that there is definitely a minority carrier in the
        %TL and so may not be as important as I think 
        %Could depend on which side of the interface the recombination
        %happens on?
        elseif surface == 1 && doped ~= 0.5
            par.sn(2) = params{i,j}(2);
            par.sn(4) = params{i,j}(2);
            par.sp(2) = params{i,j}(2);
            par.sp(4) = params{i,j}(2);
        elseif surface == 1 && doped == 0.5
            par.sn(4) = params{i,j}(2);
            par.sp(4) = params{i,j}(2);
        end

        %ion conc
        if i ~= n_ion_concs
            par.Ncat(:) = params{i,j}(1);
            par.Nani(:) = params{i,j}(1);
        end 
        
        %Do this as it seesm to reduce the discrepency between surface
        %volumetric surace recombination model and the abrupt interface one        
        %par.frac_vsr_zone = 0.05;
        par = refresh_device(par);

        soleq{i,j} = equilibrate(par);
        
        %electron only scan
        if i == n_ion_concs 
            Voc_max = 1.2;
            num_points = 281;
            while Voc_max >= Voc_max_lim
                try            
                    solCV{i, j} = doCV(soleq{i, j}.el, illumination, -0.2, Voc_max, -0.2, 1, 1, num_points);           
                    error_log(i,j) = 0;
                    results{i,j} = CVstats(solCV{i, j});
                    Voc_max = 0;                
                catch
                    if Voc_max > Voc_max_lim
                        warning("Electronic-only JV solution failed, reducing Vmax by 0.03 V")
                        Voc_max = Voc_max - 0.03;
                        num_points = num_points - 6;
                    elseif Voc_max == Voc_max_lim
                        warning("Electronic-only JV solution failed.")
                        error_log(i,j) = 1;
                        results{i,j} = 0;
                    end
                end
            end
        
        else
            Voc_max = 1.2;
            num_points = 281; 
            while Voc_max >= Voc_max_lim
                try
                    solCV{i, j} = doCV(soleq{i, j}.ion, illumination, -0.2, Voc_max, -0.2, 1e-4, 1, num_points);
                    error_log(i,j) = 0;
                    results{i,j} = CVstats(solCV{i, j});
                    Voc_max = 0;
                catch
                    if Voc_max > Voc_max_lim
                        warning("Ionic JV solution failed, reducing Vmax by 0.03 V")
                        Voc_max = Voc_max - 0.03;
                        num_points = num_points - 6;
                    elseif Voc_max == Voc_max_lim
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
Stats_array = zeros(n_ion_concs, n_recom, 5);
J_srh_result = cell(2, n_recom);
J_vsr_result = cell(2,n_recom);
V_plot = cell(2, n_recom);

for i = 1:n_ion_concs
    for j = 1:n_recom
        try
            Stats_array(i,j,1) = 1e3*results{i,j}.Jsc_f;
            Stats_array(i,j,2) = results{i,j}.Voc_f;
            Stats_array(i,j,3) = results{i,j}.FF_f;
            Stats_array(i,j,4) = results{i,j}.efficiency_f;
            x = solCV{i,j}.par.x_sub;
            loss_currents = dfana.calcr(solCV{i,j},'sub');
            Vapp = dfana.calcVapp(solCV{i,j});
            end_value = (length(Vapp)-1)/2;
            J_srh = e*trapz(x, loss_currents.srh, 2)';
            J_vsr = e*trapz(x, loss_currents.vsr, 2)';
            if i == 7
                V_plot{1,j} = Vapp;
                J_srh_result{1,j} = J_srh;
                J_vsr_result{1,j} = J_vsr;
            elseif i == 8
                V_plot{2,j} = Vapp;
                J_srh_result{2,j} = J_srh;
                J_vsr_result{2,j} = J_vsr;
            end
            J_srh_OC = interp1(Vapp(1:end_value), J_srh(1:end_value), Stats_array(i,j,2));
            J_vsr_OC = interp1(Vapp(1:end_value), J_vsr(1:end_value), Stats_array(i,j,2));
            if J_srh_OC > J_vsr_OC
                Stats_array(i,j,5) = 0;
            elseif J_srh_OC < J_vsr_OC
                Stats_array(i,j,5) = 1;
            end
        catch
            warning('No Stats')
            Stats_array(i,j,:) = 0;
        end
    end
end 

%%
figure('Name', 'JV Parameter vs Recombination vs Ion Conc', 'Position', [50 50 800 800])
Colours = parula(n_ion_concs-1);
num = 2;
labels = ["J_{SC} (mA cm^{-2})", "V_{OC} (V)", "FF", "PCE (%)"];
LegendLoc = ["northeast", "southwest", "southeast", "southeast"];
if doped == 0
    lims = [[-24 -15]; [0.8 1.2]; [0.5, 0.9]; [10 23]];
elseif doped == 1
    lims = [[-24 -15]; [0.8 1.2]; [0.5, 0.9]; [10 23]];
elseif doped == 0.5
    lims = [[-24 -15]; [0.85 1.1]; [0.5, 0.9]; [10 23]];
end
box on 
for i = 1:n_ion_concs
    if surface == 0
        if i == n_ion_concs
            semilogx(1e9*tau_SRH, Stats_array(n_ion_concs,:,5).*Stats_array(n_ion_concs,:,num), 'marker', 'x', 'Color', 'black', 'LineStyle', 'none', 'MarkerSize', 10, 'HandleVisibility', 'Off')
            semilogx(1e9*tau_SRH, (1-Stats_array(n_ion_concs,:,5)).*Stats_array(n_ion_concs,:,num), 'marker', 'o', 'Color', 'black', 'LineStyle', 'none', 'MarkerSize', 10, 'HandleVisibility', 'Off')
            semilogx(1e9*tau_SRH, Stats_array(n_ion_concs,:,num), 'marker', 'none', 'Color', 'black')
        else  
            %semilogx(1e9*tau_SRH, Stats_array(i-1,:,5).*Stats_array(i-1,:,num), 'marker', 'x', 'Color', Colours(i-1,:), 'LineStyle', 'none', 'MarkerSize', 10, 'HandleVisibility', 'Off')
            %semilogx(1e9*tau_SRH, (1-Stats_array(i-1,:,5)).*Stats_array(i-1,:,num), 'marker', 'o', 'Color', Colours(i-1,:), 'LineStyle', 'none', 'MarkerSize', 10, 'HandleVisibility', 'Off')
            semilogx(1e9*tau_SRH, Stats_array(i,:,num), 'marker', 'none', 'Color', Colours(i,:))
            hold on
        end
    elseif surface == 1
        if i == n_ion_concs
            semilogx(v_sr, Stats_array(n_ion_concs,:,5).*Stats_array(n_ion_concs,:,num), 'marker', 'x', 'Color', 'black', 'LineStyle', 'none', 'MarkerSize', 10, 'HandleVisibility', 'Off')
            semilogx(v_sr, (1-Stats_array(n_ion_concs,:,5)).*Stats_array(n_ion_concs,:,num), 'marker', 'o', 'Color', 'black', 'LineStyle', 'none', 'MarkerSize', 10, 'HandleVisibility', 'Off')
            semilogx(v_sr, Stats_array(n_ion_concs,:,num), 'marker', 'none', 'Color', 'black')
        else  
            %semilogx(v_sr, Stats_array(i-1,:,5).*Stats_array(i-1,:,num), 'marker', 'x', 'Color', Colours(i-1,:), 'LineStyle', 'none', 'MarkerSize', 10, 'HandleVisibility', 'Off')
            %semilogx(v_sr, (1-Stats_array(i-1,:,5)).*Stats_array(i-1,:,num), 'marker', 'o', 'Color', Colours(i-1,:), 'LineStyle', 'none', 'MarkerSize', 10, 'HandleVisibility', 'Off')
            semilogx(v_sr, Stats_array(i,:,num), 'marker', 'none', 'Color', Colours(i,:))
            hold on
        end
    end
end
set(gca, 'Fontsize', 25)
if surface == 0
    xlabel('Shockley-Read-Hall lifetime (ns)', 'FontSize', 30)
    xlim([1, 1000])
elseif surface == 1
    xlabel('Surface Recombination Velocity (cm s^{-1})', 'FontSize', 30)
    xlim([1e-1, 1e2])
end
ylabel(labels(num), 'FontSize', 30)
ylim(lims(num,:))
legend({' 5e15', ' 1e16', ' 5e16', ' 1e17', ' 5e17', ' 1e18', ' No Ions'}, 'Location', 'southeast', 'FontSize', 25, 'NumColumns', 2)
title(legend, 'Ion Concentration (cm^{-3})', 'FontSize', 25)

%% Plot JV as a function of recombination parameter for a given ion conc
figure('Name', 'JVPlot', 'Position', [100 100 800 800])
Colours = parula(n_recom);
%Ion_Conc = [1e15 5e15 1e16 5e16 1e17 5e17 1e18 0];
%Set which ion concentration to plot for
%Have coppied the array above so you can see which nuber is the right one
%easily
num_Ion_Conc = 6;

for j = 1:n_recom
    v = dfana.calcVapp(solCV{num_Ion_Conc, j});
    J = dfana.calcJ(solCV{num_Ion_Conc, j}).tot(:,1);
    
    hold on
    xline(0, 'black', 'HandleVisibility', 'off')
    yline(0, 'black', 'HandleVisibility', 'off')
    plot(v(:), J(:)*1000, 'color', Colours(j,:), 'LineWidth', 3) 
    hold off

end

box on 
set(gca, 'FontSize', 25)
xlim([-0.15, 1.2])
ylim([-25,5])
if surface == 0
    legend({'10', '30', '50', '100', '300', '500', '1000'}, 'Location', 'northwest', 'FontSize', 25, 'NumColumns', 2)
    title(legend, 'SRH Lifetime (ns)', 'FontSize', 25)
elseif surface == 1
    legend({'0.1', '1', '10', '100', '1000', '10000', '100000'}, 'Location', 'eastoutside', 'FontSize', 25, 'NumColumns', 2)
    title(legend, 'Surface Recombination\newlineVelocity (cm s^{-1})', 'FontSize', 25)
end
xlabel('Voltage(V)', 'FontSize', 30)
ylabel('Current Density (mAcm^{-2})', 'FontSize', 30)
ax1 = gcf;

%% Save results and solutions
save_file = 1;
if save_file == 1
    if doped == 0 && surface == 0
        filename = 'v5_undoped_tauSRH.mat';
    elseif doped == 0 && surface == 1
        filename = 'v5_undoped_vsr.mat';
    elseif doped == 1 && surface == 0
        filename = 'v5_doped_tauSRH.mat';
    elseif doped == 1 && surface == 1
        filename = 'v5_doped_vsr.mat';
    elseif doped == 0.5 && surface == 1
        filename = 'SaP_params_vsr.mat';
    end 
    save(filename, 'results', 'solCV')
end
