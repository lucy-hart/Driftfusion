par = pc('Input_files/SAM_MAFACsPbIBr_C60.csv');
par.z_a = 0;
par.vsr_mode = 0;
par = refresh_device(par);
eqm = equilibrate(par);
sol_light = lightonRs(eqm.ion,1,1,1,0,100);
%JVsol = doCV(eqm.el, 1, -0.2, 1.2, -0.2, 10e-3, 1, 281);