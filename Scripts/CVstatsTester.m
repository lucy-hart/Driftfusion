par = pc('Input_files/pedotpss_mapi_pcbm.csv');

%sol = equilibrate(par);

sol_CV1 = doCV(sol.ion, 1, 1.0, 1.1, 0, 100e-3, 3, 241);
sol_CV2 = doCV(sol.ion, 1, 0, 1.2, 0, 100e-3, 3, 241);
sol_CV3 = doCV(sol.ion, 1, 0, 1.2, 0, 100e-3, 1, 241);