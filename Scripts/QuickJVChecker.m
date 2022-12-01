%par=pc('Input_files/EnergyOffsetSweepParameters_v2.csv');
%par = pc('Input_files/PTAA_MAPI_NegOffset.csv');
par = pc('Input_files/PTAA_MAPI_NoOffset.csv');
%par = pc('Input_files/PTAA_MAPI_PosOffset.csv');

Fiddle_with_Energetics = 0;
%%
if Fiddle_with_Energetics == 1

    %row
    DHOMO = 0.3;
    %DHOMO = Delta_HOMO(4);
    %column
    DLUMO = -0.0;
    %DLUMO = Delta_LUMO(11);

    %HTL Energetics
    par.Phi_left = -5.3;
    par.Phi_IP(1) = par.Phi_IP(3) + DHOMO;
    par.Phi_EA(1) = par.Phi_IP(1) + 2.5;
    par.EF0(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
    par.Et(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
    if par.Phi_left < par.Phi_IP(1)
        par.Phi_left = par.Phi_IP(1);
    end
    
    %ETL Energetics
    par.Phi_right = -4.1;
    par.Phi_EA(5) = par.Phi_EA(3) + DLUMO;
    par.Phi_IP(5) = par.Phi_EA(5) - 2.5;
    par.EF0(5) = (par.Phi_IP(5)+par.Phi_EA(5))/2;
    par.Et(5) = (par.Phi_IP(5)+par.Phi_EA(5))/2;
    if par.Phi_right > par.Phi_EA(5)
        par.Phi_right = par.Phi_EA(5);
    end

    par = refresh_device(par);

end

eqm_QJV = equilibrate(par);

%%
CV_sol_ion = doCV(eqm_QJV.ion, 1, -0.2, 1.25, -0.2, 1e-4, 1, 291);
%CV_sol_el = doCV(eqm_QJV.el, 1, -0.2, 1.17, -0.2, 1e-4, 1, 275);

%%
Plot_Current_Contributions(CV_sol_ion)
stats_ion = CVstats(CV_sol_ion)
%stats_el = CVstats(CV_sol_el)


