num = 1;

OneSunSol = changeLight(eqm_solutions_dark{num}.ion, 1, 0);
VappStruct = genVappStructs(OneSunSol, [0.6, Voc_ion(num), 1.2], 1);