par = pc('Input_files/SingleLayer_SurfaceTraps.csv');
par.vsr_mode = 0;

par = refresh_device(par);
eqm = equilibrate(par,0,1);
eqm.ion.par.AbsTol = 1e-3;
eqm.ion.par.RelTol = 1e-6;
%% 
eqm.ion.par.Rs = 1e6;
sollight = changeLight(eqm.ion,1,1e-3);
%%
sollight.par.tmax = 1;
sollight.par.tmesh_type = 'linear';
sol_light = df(sollight);