
par=pc('Input_files/EnergyOffsetSweepParameters_v5_undoped_Weidong_fiddled.csv');
par.light_source1 = 'laser';
par.laser_lambda1 = 532;
par.pulsepow = 62;
par.RelTol_vsr = 0.1;
par = refresh_device(par);
eqm_QJV = equilibrate(par);
%%
suns = 1;
JV_sol_ion = doCV(eqm_QJV.ion, suns, -0.2, 1.2, -0.2, 10e-3, 1, 281);
Plot_Current_Contributions(JV_sol_ion)
stats_ion = CVstats(JV_sol_ion)

