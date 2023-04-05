%par=pc('Input_files/EnergyOffsetSweepParameters_v2.csv');
%par=pc('Input_files/EnergyOffsetSweepParameters_v3.csv');
par = pc('Input_files/PTAA_MAPI_NegOffset.csv');
%par = pc('Input_files/PTAA_MAPI_NegOffset_lowerVbi.csv');
%par = pc('Input_files/PTAA_MAPI_NoOffset.csv');
%par = pc('Input_files/PTAA_MAPI_PosOffset.csv');
%par = pc('Input_files/PTAA_MAPI_PCBM_ForPaper.csv');

Fiddle_with_Energetics = 0;
Fiddle_with_IonConc = 0;
IonConc = 1e15;
%%
if Fiddle_with_Energetics == 1

    %row
    DHOMO = -0.15;
    %DHOMO = Delta_HOMO(4);
    %column
    DLUMO = 0.15;
    %DLUMO = Delta_LUMO(11);

    %HTL Energetics
    par.Phi_left = -5.15;
    par.Phi_IP(1) = par.Phi_IP(3) + DHOMO;
    par.Phi_EA(1) = par.Phi_IP(1) + 2.5;
    par.EF0(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
    par.Et(1) = (par.Phi_IP(1)+par.Phi_EA(1))/2;
    if par.Phi_left < par.Phi_IP(1)
        par.Phi_left = par.Phi_IP(1);
    end
    
    %ETL Energetics
    par.Phi_right = -4.05;
    par.Phi_EA(5) = par.Phi_EA(3) + DLUMO;
    par.Phi_IP(5) = par.Phi_EA(5) - 2.5;
    par.EF0(5) = (par.Phi_IP(5)+par.Phi_EA(5))/2;
    par.Et(5) = (par.Phi_IP(5)+par.Phi_EA(5))/2;
    if par.Phi_right > par.Phi_EA(5)
        par.Phi_right = par.Phi_EA(5);
    end

    par = refresh_device(par);

end

if Fiddle_with_IonConc == 1

   par.Ncat(:) = IonConc;
   par.Nani(:) = IonConc;

   par = refresh_device(par);

end

eqm_QJV = equilibrate(par);

%%
CV_sol_ion = doCV(eqm_QJV.ion, 1.1, -0.2, 1.2, -0.2, 1e-4, 1, 281);
%CV_sol_el = doCV(eqm_QJV.el, 1.1, -0.2, 1.25, -0.2, 1e-4, 1, 291);

Plot_Current_Contributions(CV_sol_ion) 
stats_ion = CVstats(CV_sol_ion)
%stats_el = CVstats(CV_sol_el)

%%
%Make one sun solution at a given applied voltage
run = 0; 
if run == 1
    %Vapp = 0.62; %Uniform ion distribution from JV for negtive offset case
    Vapp = stats_ion.Voc_f; %Uniform ion distribution from JV
    sol_ill = changeLight(eqm_QJV.ion, 1.1, 0, 1);
    sol_ill_bias = genVappStructs(sol_ill, Vapp, 1);
    try
        CV_sol_startbias = doCV(sol_ill_bias, 1.1, Vapp, -0.20, 1.20, 100, 1, 281);
    catch
        warning('No joy.')
    end

    Plot_Current_Contributions(CV_sol_startbias)
    stats_ions_bias = CVstats(CV_sol_startbias)

    dfplot.ELxnpxacx(CV_sol_ion, 1e4*(0.2+Vapp))
    dfplot.ELxnpxacx(sol_ill_bias, sol_ill_bias.t(end))
end     


