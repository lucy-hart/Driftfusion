%% m-i-m for Noushin 

par=pc('Input_files/1_layer_test.csv');

% Work Functions: Au = -5.1 eV, Ag = -4.7 eV, Al = -4.3 eV, ITO = many and
% varied, found -4.3 eV, -4.6 eV, 4.2 eV and -4.8 eV

anode_WF = -4.7;
cathode_WF = -4.3;

par.Phi_left = anode_WF;
par.Phi_right = cathode_WF;

eqm_QJV = equilibrate(par);

%%
suns = 0;
V_max = 1.2;
V_min = -0.2;
scan_rate = 0.01;

JV_sol = doCV(eqm_QJV.ion, suns, V_min, V_max, V_min, scan_rate, 1, 200*(V_max-V_min)+1);

if suns == 1
    Plot_Current_Contributions_v2(JV_sol) 
    stats_ion = CVstats(JV_sol)
end 



