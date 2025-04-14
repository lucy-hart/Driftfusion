par=pc('Input_files/EnergyOffsetSweepParameters_v5_undoped_Weidong.csv');

ions = [1e15 1e16 1e17 1e18];
eqms = cell(1,length(ions));
ill =  cell(1,length(ions));
for i = 1:length(ions)
    par.Nani(:) = ions(i); 
    par.Ncat(:) = ions(i);
    par = refresh_device(par);
    eqms{i} = equilibrate(par);
    ill{i} = changeLight(eqms{i}.ion, 1, 0, 1);
end