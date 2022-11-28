par=pc('Input_files/EnergyOffsetSweepParameters.csv');
%par = pc('Input_files/PTAA_MAPI_NegOffset.csv');
%par = pc('Input_files/PTAA_MAPI_NoOffset.csv');
%par = pc('Input_files/PTAA_MAPI_PosOffset.csv');

DHOMO = 0.3;
DLUMO = -0;

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

eqm_QJV = equilibrate(par);

%%
CV_sol_ion = doCV(eqm_QJV.ion, 1, -0.2, 1.15, -0.2, 1e-4, 1, 271);
CV_sol_el = doCV(eqm_QJV.el, 1, -0.2, 1.18, -0.2, 10e-3, 1, 277);

%%
Plot_Current_Contributions(CV_sol_ion)
stats_ion = CVstats(CV_sol_ion)
stats_el = CVstats(CV_sol_el)


