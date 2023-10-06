data = readtable('C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Davide_OrganicPeroHybrid\EQE_data_minus0p5V.xlsx');
%%
h = 6.63e-34;
c = 3e8;
k_B = 1.38e-23;
T = 298;
q = 1.6e-19;

wl_stop = 1070;
num_points = (wl_stop-300)/10 + 1;

lambda = linspace(300e-9, wl_stop*1e-9, num_points);
E = h*c./(lambda);
E_in_eV = E/q;

wl = data{2:end,1};
[~,argstop] = min(abs(wl-wl_stop));
PM6 = data{2:argstop+1,3}/100;
PM7 = data{2:argstop+1,4}/100;
PCE12 = data{2:argstop+1,5}/100;

phi_E = (2*pi/(h^3*c^2)).*(E.^2./(exp((E)/(k_B*T))-1));

%Factor of two as two side to the solar cell 
%Units are A cm-2 (I think - that's why I have factor of 10^4)
J_rad_PM6 = 2*1e4*q*trapz(flip(E), flip(PM6.*phi_E'));
J_rad_PM7 = 2*1e4*q*trapz(flip(E), flip(PM7.*phi_E'));
J_rad_PCE12 = 2*1e4*q*trapz(flip(E), flip(PCE12.*phi_E'));

%%
figure('Name', 'log EQE plot')
semilogy(flip(E/1.6e-19), flip(PM6))
hold on 
semilogy(flip(E/1.6e-19), flip(PM7))
semilogy(flip(E/1.6e-19), flip(PCE12))

xlim([1, 2.5])
ylim([1e-5, 1])
legend({'PM6:Y6', 'PM7:Y6', 'PCE12:Y6'}, 'Location', 'southeast')
