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
%
symmetric = 1;
mu = 1;
Nion = 1e17;
eps_pero = 25;

%
mu_TL_ar = logspace(log10(5e-6),-3,7);
n_TL_mus = length(mu_TL_ar);

if mu == 1
    mu_pero = logspace(log10(0.1),log10(50),7);
    n_pero_var  = length(mu_pero);
    %Rows are the TL mobilities    
    %Columns are the perovskite mobilities
    params = cell(n_TL_mus, n_pero_var);
    for i=1:n_TL_mus
        for j=1:n_pero_var
            params{i,j} = [mu_TL_ar(i), mu_pero(j)];
        end
    end
elseif mu == 0
    tau_SRH = [1e-8 3e-8 5e-8 1e-7 3e-7 5e-7 1e-6];
    n_pero_var  = length(tau_SRH);
    %Rows are the TL mobilities    
    %Columns are the SRH lifetimes 
    params = cell(n_TL_mus, n_pero_var);
    for i=1:n_TL_mus
        for j=1:n_pero_var
            params{i,j} = [mu_TL_ar(i), tau_SRH(j)];
        end
    end
end

%%
%Select the correct input file for doped or undoped cases
if symmetric == 0
    par=pc('Input_files/EnergyOffsetSweepParameters_v5_undoped_Weidong.csv');
    Voc_max_lim = 1.05;
elseif symmetric == 1
    par=pc('Input_files/EnergyOffsetSweepParameters_v5_undoped.csv');
    Voc_max_lim = 1.05;
end

%% Choose the energetics of the TLs 
%Default values are FILL THIS IN
%Will use these values if Fiddle_with_Energetics is 0
Fiddle_with_Energetics = 1;

if Fiddle_with_Energetics == 1
    %Choose the offsets for the system
    %Positive offset for DHOMO means TL VB lies above the perovskite VB
    %Negative offset for DLUMO means TL CB lies below the perovskite CB
    DHOMO = 0;%.2;
    DLUMO = -0;%.2;
    if symmetric == 1
        %HTL Energetics
        par.Phi_left = -5.15;
        par.Phi_IP(1) = par.Phi_IP(3) + DHOMO;
        par.Phi_EA(1) = par.Phi_IP(1) + 2.5;
        par.EF0(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
        par.Et(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
        if par.Phi_left < par.Phi_IP(1) + 0.1
            par.Phi_left = par.Phi_IP(1) + 0.1;
        end
        %ETL Energetics
        par.Phi_right = -4.05;
        par.Phi_EA(5) = par.Phi_EA(3) + DLUMO;
        par.Phi_IP(5) = par.Phi_EA(5) - 2.5;
        par.EF0(5) = (par.Phi_IP(5)+par.Phi_EA(5))/2;
        par.Et(5) = (par.Phi_IP(5)+par.Phi_EA(5))/2;
        if par.Phi_right > par.Phi_EA(5) - 0.1
            par.Phi_right = par.Phi_EA(5) - 0.1;
        end
    elseif symmetric == 0
        %ETL Energetics
        par.Phi_right = -4.00;
        par.Phi_EA(5) = par.Phi_EA(3) + DLUMO;
        par.Phi_IP(5) = par.Phi_EA(5) - 2.5;
        par.Et(5) = (par.Phi_IP(5) + par.Phi_EA(5))/2;
        par.EF0(5) = (par.Phi_IP(5) + par.Phi_EA(5))/2;
        if par.Phi_right > par.Phi_EA(5) - 0.05
            par.Phi_right = par.Phi_EA(5) - 0.05;
        end
    end
    par.Ncat(:) = Nion;
    par.Nani(:) = Nion;
    par.epp(2:4) = eps_pero;
    par = refresh_device(par);

end

%% Set up structures for storing the results
error_log = zeros(n_TL_mus, n_pero_var);
soleq = cell(n_TL_mus, n_pero_var);
solCV = cell(n_TL_mus, n_pero_var);
results = cell(n_TL_mus, n_pero_var);

%% Do (many) JV sweeps

%Set the illumnation for the JV sweeps 
illumination = 1;

%Reset the electrode work functions in each loop to be safe as they are
%changed for the cases where E_LUMO (E_HOMO) is far below (above) the CB
%(VB)
for i = 1:n_TL_mus
    for j = 1:n_pero_var
        if mu == 0
            disp(["tau_SRH = ", num2str(tau_SRH(j)), " s"])
        elseif mu == 1
            disp(["mu_pero = ", num2str(mu_pero(j)), " cm2 V-1 s-1"])
        end 

        if mu == 0
            par.taun(3) = params{i,j}(2);
            par.taup(3) = params{i,j}(2);

        elseif mu == 1 
            par.mu_n(3) = params{i,j}(2);
            par.mu_p(3) = params{i,j}(2);
        end

        if symmetric == 1
            par.mu_n(1) = params{i,j}(1);
            par.mu_p(1) = params{i,j}(1);
        end
        par.mu_n(5) = params{i,j}(1);
        par.mu_p(5) = params{i,j}(1);

        par.light_source1 = 'laser';
        par.laser_lambda1 = 532;
        par.pulsepow = 62;
        par.RelTol_vsr = 0.1;
        par = refresh_device(par);

        soleq{i,j} = equilibrate(par);
       
        Voc_max = 1.3;
        num_points = 301; 
        while Voc_max >= Voc_max_lim
            try
                solCV{i, j} = doCV(soleq{i, j}.ion, illumination, -0.2, Voc_max, -0.2, 10e-3, 1, num_points);
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

toc

%% Extract results 
Stats_array = zeros(n_TL_mus, n_pero_var, 5);

num_start = sum(solCV{1,6}.par.layer_points(1:2))+1;
num_stop = num_start + solCV{1,6}.par.layer_points(3)-1;
x = solCV{1,6}.par.x_sub;
d = solCV{1,6}.par.d(3);

for i = 1:n_TL_mus
    for j = 1:n_pero_var
        try
            Stats_array(i,j,1) = 1e3*results{i,j}.Jsc_f;
            Stats_array(i,j,2) = results{i,j}.Voc_f;
            Stats_array(i,j,3) = results{i,j}.FF_f;
            Stats_array(i,j,4) = results{i,j}.efficiency_f;
            %Calculate QFLS 
            [~, ~, Efn_ion, Efp_ion] = dfana.calcEnergies(solCV{i,j});
            QFLS_ion = trapz(x(num_start:num_stop), Efn_ion(:, num_start:num_stop)-Efp_ion(:,num_start:num_stop),2)/d;            
            Stats_array(i,j,5) = QFLS_ion(21,:);

        catch
            warning('No Stats')
            Stats_array(i,j,:) = 0;
        end
    end
end 

%%
figure('Name', 'JV Parameter vs Recombination vs Ion Conc', 'Position', [50 50 800 800])
Colours = parula(n_TL_mus);
num = 5;
labels = ["J_{SC} (mA cm^{-2})", "V_{OC} (V)", "FF", "PCE (%)", "QFLS_{SC}"];
LegendLoc = ["northeast", "southwest", "southeast", "southeast"];
if symmetric == 0
    lims = [[-24 -5]; [0.8 1.2]; [0.5, 0.9]; [10 23]; [0.8 1.1]];
elseif symmetric == 1
    lims = [[-24 -5]; [0.8 1.2]; [0.5, 0.9]; [10 23]; [0.8 1.1]];
elseif symmetric == 0.5
    lims = [[-24 -5]; [0.85 1.1]; [0.5, 0.9]; [10 23]; [0.8 1.1]];
end
box on 

for i = 1:n_TL_mus
    if mu == 0
       semilogx(1e9*tau_SRH, Stats_array(i,:,num))
       hold on
    elseif mu == 1
        semilogx(mu_pero, Stats_array(i,:,num))
        hold on
    end
end
hold off
set(gca, 'Fontsize', 25)
if surface == 0
    xlabel('Shockley-Read-Hall lifetime (ns)', 'FontSize', 30)
    xlim([1, 1000])
elseif surface == 1
    xlabel('Perovskite Mobility (cm^{2} V^{-1} s^{-1})', 'FontSize', 30)
    xlim([1e-1, 50])
end
ylabel(labels(num), 'FontSize', 30)
ylim(lims(num,:))
legend({' 1e-3', ' 1e-4', ' 1e-5', ' 1e-6'}, 'Location', 'southeast', 'FontSize', 25)
title(legend, 'Transport Layer Mobility (cm^{2} V^{-1} s^{-1})', 'FontSize', 25)

%% Save results and solutions
save_file = 0;
if save_file == 1
    filename = 'SaP_params_vsr.mat';
    save(filename, 'results', 'solCV')
end
