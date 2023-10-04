%par_test = pc('Input_files/SAM_MAPI_PM6Y6_NoSurf_CheckDarkCurrent.csv');
par_test = pc('Input_files/SAM_MAPI_C60.csv');
eqm_test = equilibrate(par_test);
voltage_ar = linspace(-0.5, 0.2, 15);
Jdark = doDarkJV(eqm_test.ion, voltage_ar, 5);