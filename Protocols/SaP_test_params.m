function Jvalues = SaP_test_params(params, Vbias, Vpulse)
    %Takes in a params variable which contains the values to test for this run
    %in the optimisation process. Should be in the order
    %params = [IonConc, EC_HTL, EV_ETL, EF_HTL, EF_ETL,
    %          tau_n, tau_p, vsurf_HTL, vsurf_ETL, 
    %          mu_n,pero, mu_p,pero, mu_n,ETL, mu_p,HTL,
    %          epsilon_ETL] 
    %Returns the J values for each Vpulse and each Vbias in the form of a 
    %[num_bias num_pulse] shaped matrix
    
    params = cell2mat(params);
    Vbias = cell2mat(Vbias);
    Vpulse = cell2mat(Vpulse);

    par=pc('Input_files/TiO2_MAPI_Spiro_TestSaP.csv');

    par.vsr_check = 0;
    
    %Mobile Ion Concentration
    par.Ncat(:) = params(1);
    par.Nani(:) = params(1);
               
    %HTL Energetics 
    par.Phi_IP(1) = par.Phi_IP(3) + params(2);
    par.Phi_EA(1) = par.Phi_IP(1) + 2.5;
    par.EF0(1) = par.Phi_IP(1) + params(4);
    par.Et(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
    par.Phi_left = par.EF0(1);
    
    %ETL Energetics
    par.Phi_EA(5) = par.Phi_EA(3) - params(3);
    par.Phi_IP(5) = par.Phi_EA(5) - 2.5;
    par.EF0(5) = par.Phi_EA(5) - params(5);
    par.Et(5) = (par.Phi_IP(5) + par.Phi_EA(5))/2;
    par.Phi_right = par.EF0(5);
    
    %Bulk recombination 
    par.taun(2:4) = params(6);
    par.taup(2:4) = params(7);
    
    %Surface Recombination (HTL)
    par.sn(2) = params(8);
    par.sp(2) = params(8);
    
    %Surface Recombination (ETL)
    par.sn(4) = params(9);
    par.sn(4) = params(9);
    
    %Mobilities
    par.mu_n(3) = params(10);
    par.mu_p(3) = params(11);
    par.mu_n(5) = params(12);
    par.mu_p(1) = params(13);
    
    %Epsilon (ETL only in case this matters for SnO2 versus TiO2)
    par.epp(5) = params(14);
    
    try
        eqm = equilibrate(par);
    catch
        warning('Could not equilibrate device for these input parameters.')
        eqm = 0;
    end
    
    %% Do the SaP measurement
    num_bias = length(Vbias);
    num_pulse = length(Vpulse);

    if not(isa(eqm, 'struct'))
        Jvalues = zeros([num_bias num_pulse]);

    else 
        suns = 1;
        tramp = 8e-4;
        tsample = 1e-3;
        tstab = 200;
        
        sol = doSaP_multi(eqm.ion, Vbias, Vpulse, tramp, tsample, tstab, suns, 0);
        
        %% Extract the current values 
        Jvalues = zeros([num_bias num_pulse]);
        for i = 1:num_bias
            for j = 1:num_pulse
                Jvalues(i,j) = sol{i,j+1}.Jpulse;
            end
        end
    end
end