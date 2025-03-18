%% Define parameter space
%Choose to use doped or undoped TLs
doped = 0;

n_values = 53;
params = zeros(n_values,5);

Delta_TL = 0.4;
Ion_Conc = 3;
mu_pero = 3;
mu_ETL = 3;
%Lamb = [800e-7 1200e-7 2000e-7];
tau_pero = 3;

%seeds: 0, 42, 666
rng(666,'twister');

% idxes = randi(length(Lamb),n_values,1);
params(:,1) = -0.1 + Delta_TL.*rand(n_values,1);
params(:,2) = 10.^(15 + Ion_Conc.*rand(n_values,1));
params(:,3) = 10.^(-1 + mu_pero.*rand(n_values,1));
params(:,4) = 10.^(-5 + mu_ETL.*rand(n_values,1));
params(:,5) = 10.^(-9 + tau_pero.*rand(n_values,1));
% l_amb = Lamb(idxes)';
% params(:,5) = 0.5*(l_amb.^2)./(0.0257.*params(:,3));

% Set up structures for storing the results
results = zeros(n_values,3);

%% Do (many) JV sweeps
%Select the correct input file for doped or undoped cases 
if doped == 1
    par=pc('Input_files/EnergyOffsetSweepParameters_v5_doped.csv');
elseif doped == 0
    par=pc('Input_files/EnergyOffsetSweepParameters_v5_undoped_Weidong.csv');
end 

%Set the illumination for the JV sweeps and the max voltage to go up to in
%the JV sweeps

illumination = 1;
max_val = 1.2;
phi_R = -4.0;

%Reset the electrode work functions in each loop to be safe as they are
%changed for the cases where E_LUMO (E_HOMO) is far below (above) the CB
%(VB)
for i = 1:n_values
    disp(num2str(i))

    %ETL Energetics
    %Need to use opposite sign at ETL to keep energy offsets symmetric
    par.Phi_right = phi_R;
    par.Phi_EA(5) = par.Phi_EA(3) - params(i,1);
    par.Phi_IP(5) = par.Phi_EA(5) - 2.5;
    par.Et(5) = (par.Phi_IP(5) + par.Phi_EA(5))/2;
    par.EF0(5) = (par.Phi_IP(5) + par.Phi_EA(5))/2;
    if par.Phi_right > par.Phi_EA(5) - 0.05
        par.Phi_right = par.Phi_EA(5) - 0.05;
    end

    %ion conc

    par.Ncat(:) = params(i,2);
    par.Nani(:) = params(i,2);

    par.mu_n(3) = params(i,3);
    par.mu_p(3) = params(i,3);
    par.mu_n(5) = params(i,4);
    par.mu_p(5) = params(i,4);

    par.taun(3) = params(i,5);
    par.taup(3) = params(i,5);

    par.light_source1 = 'laser';
    par.laser_lambda1 = 532;
    par.pulsepow = 62;
    par.RelTol_vsr = 0.1;
    par = refresh_device(par);

    soleq  = equilibrate(par);
    
    Fermi_offset = par.EF0(5) - par.EF0(1);
    if Fermi_offset > max_val
        Voc_max = Fermi_offset + 0.05;
    else
        Voc_max = max_val;
    end
    num_points = (2*100*(Voc_max+0.2))+1; 
    while Voc_max >= 1.05
        try
            solCV = doCV(soleq.ion, illumination, -0.2, Voc_max, -0.2, 10e-3, 1, num_points);
            num_start = sum(solCV.par.layer_points(1:2))+1;
            num_stop = num_start + solCV.par.layer_points(3)-1;
            x = solCV.par.x_sub;
            d = solCV.par.d(3);
            [~, ~, Efn_ion, Efp_ion] = dfana.calcEnergies(solCV);
            QFLS_ion = trapz(x(num_start:num_stop), Efn_ion(:, num_start:num_stop)-Efp_ion(:,num_start:num_stop),2)/d;            
            results(i,1) = QFLS_ion(21,:);
            stats = CVstats(solCV);
            results(i,2) = stats.Jsc_f;
            results(i,3) = stats.FF_f;
            Voc_max = 0;
        catch
            if Voc_max > 1.05
                warning("Ionic JV solution failed, reducing Vmax by 0.03 V")
                Voc_max = Voc_max - 0.03;
                num_points = num_points - 6;
            elseif Voc_max == 1.05
                warning("Ionic JV solution failed.")
                results(i,:) = 0;
            end
        end
    end
end

