par=pc('Input_files/HighEfficiencyPaper-ForBandDiagrams.csv');
%par.z_c = -1;

Ef_pero = linspace(-4,-4.6, 7);
%Ef_C60 = linspace(-4.3,-4.9,7);
%Ef_C60 = linspace(-4.1,-4.7,7);
Ec_plot = zeros(1101, length(Ef_pero));
Ev_plot = zeros(1101, length(Ef_pero));
for i = 1:length(Ef_pero)
    par.EF0(1:2) = Ef_pero(i);
    %par.EF0(3) = Ef_C60(i);
    par = refresh_device(par);
    eqm = equilibrate(par);
    [Ec, Ev, Efn, Efp] = dfana.calcEnergies(eqm.ion);
    Ec_plot(:,i) = Ec(200,:)';
    Ev_plot(:,i) = Ev(200,:)';
end
x = eqm.ion.x';