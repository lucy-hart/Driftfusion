function plot_tail_state_schematic(Ef, E_UCB, Nt0, Nc, Nv, n, p, Cn_CB, Cp_CB)

kB = 8.617330350e-5;
T = 300;
E = -1:0.01:0;
Ecb = E(end);
Evb = E(1);
E_discrete = -1:0.1:0;
FD_prob = 1./(exp((E - Ef)./(kB*T)) + 1);

Nt = Nt0*exp((E - Ef)./E_UCB);
nt = Nt.*(Cn_CB*n + Cp_CB*Nv.*exp((Evb - E)./(kB*T)))./...
    (Cn_CB.*(n + Nc.*exp((E - Ecb)./(kB*T))) + Cp_CB.*(p + Nv.*exp((Evb - E)./(kB*T))));

nt_discrete = Nt0*exp((E_discrete - Ef)./E_UCB).*(Cn_CB*n + Cp_CB*Nv.*exp((Evb - E_discrete)./(kB*T)))./...
    (Cn_CB.*(n + Nc.*exp((E_discrete - Ecb)./(kB*T))) + Cp_CB.*(p + Nv.*exp((Evb - E_discrete)./(kB*T))));

figure(700)
yyaxis right
area(E, nt)
hold on
semilogy(E, Nt, E_discrete, nt_discrete)
hold off
set(gca, 'YScale', 'log')
ylabel('Density per unit energy (cm-3 eV-1)')
xlabel('Energy (eV)')
yyaxis left
plot(E, FD_prob, 'k--')
ylabel('Occupation probability')
end