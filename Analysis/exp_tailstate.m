function exp_tailstate(par, n, p, N_levels)
% Written for a single layer device only

kT = par.kB*par.T;
Ecb = 0;
Evb = -1.6;
E_UCB = 3*kT;         % Characteristic energy;
Nt_CBedge = 1e20;       % Density of trap states at CB edge
E = Evb:(Ecb - Evb)/(N_levels-1):Ecb;

Cn_CB = par.Cn_CB(1);
Cp_CB = par.Cp_CB(1);
ni = par.ni(1);
Nc = par.Nc(1);
Nv = par.Nv(1);

Nt_CBTS_continuous = Nt_CBedge*exp(Ecb - (-E_UCB/kT));
for j = 1:N_levels
    Nt_CB(j) = Nt_CBedge*exp((E(j) - Ecb)/E_UCB);
    nt_CB(j) = Nc*exp((E(j) - Ecb)/kT);
    pt_CB(j) = Nv*exp((Evb - E(j))/kT);
    
    n_trap(j) = Nt_CB(j)*((Cn_CB*n + Cp_CB*pt_CB(j))/(Cn_CB*(n + nt_CB(j)) + Cp_CB*(p + pt_CB(j))));

    r_srh_CBTS(j) = (Nt_CB(j)*Cn_CB*Cp_CB*(n*p - ni^2))/(Cp_CB*(p + pt_CB(j)) + Cn_CB*(n + nt_CB(j)));
end

r_srh_CBTS_element = Nt_CBedge.*exp((E(j) - Ecb)./E_UCB).*(Cn_CB.*Cp_CB*(n*p - ni^2))./(Cp_CB*(p + pt_CB) + Cn_CB*(n + nt_CB));
r_srh_CBTS_sum = sum(r_srh_CBTS_element, 2)
r_srh_CBTS_int = integral(@r_srh_CBTS_fun, Evb, Ecb)
r_srh_CBTS_integral_external_fun = calc_rsrh_CBTS(Evb, Ecb, Evb, Ecb, E_UCB, Nt_CBedge, Cn_CB, Cp_CB, Nc, Nv, n, p, ni, par.kB, par.T)
r_srh_VBTS_integral_external_fun = calc_rsrh_VBTS(Evb, Ecb, Evb, Ecb, E_UCB, Nt_CBedge, Cn_CB, Cp_CB, Nc, Nv, n, p, ni, par.kB, par.T)

ntrap_CBTS_integral_external_fun = calc_ntrap_CBTS(Evb, Ecb, Evb, Ecb, E_UCB, Nt_CBedge, Cn_CB, Cp_CB, Nc, Nv, n, p, ni, par.kB, par.T)
ptrap_VBTS_integral_external_fun = calc_ptrap_VBTS(Evb, Ecb, Evb, Ecb, E_UCB, Nt_CBedge, Cn_CB, Cp_CB, Nc, Nv, n, p, ni, par.kB, par.T)

%% plots
figure(300)
semilogy(E, Nt_CB, 'o')
xlabel('Energy (eV)')
ylabel('Trap state density (cm-3)')

figure(301)
semilogy(E, n_trap, 'o')
xlabel('Energy (eV)')
ylabel('Trapped electron density (cm-3)')

figure(302)
semilogy(E, r_srh_CBTS, 'o')
xlabel('Energy (eV)')
ylabel('Recombination rate (cm-3s-1)')

    function r_srh_CBTS_integrand = r_srh_CBTS_fun(E)
        Nt_CB = Nt_CBedge*exp((E - Ecb)/E_UCB);
        nt_CB = Nc*exp((E - Ecb)/kT);
        pt_CB = Nv*exp((Evb - E)/kT);
        r_srh_CBTS_integrand = (Nt_CB.*Cn_CB.*Cp_CB.*(n.*p - ni.^2))./(Cp_CB.*(p + pt_CB) + Cn_CB.*(n + nt_CB)); 
    end
end
