function r_srh_CBTS = calc_rsrh_CBTS(E_min, E_max, Evb, Ecb, E_UCB, Nt_CB, Cn_CB, Cp_CB, Nc, Nv, n, p, ni, kB, T)
% Calculates the conduction band SRH recombination rate between limits
% E_MIN and E_MAX
r_srh_CBTS_fun = @(E) (Nt_CB.*exp((E - Ecb)./E_UCB)).*(Cn_CB.*Cp_CB.*(n*p - ni^2)./...
    (Cp_CB.*(p + Nv.*exp((Evb - E)./(kB*T))) + Cn_CB.*(n + Nc.*exp((E - Ecb)./(kB*T)))));
r_srh_CBTS = integral(r_srh_CBTS_fun, E_min, E_max);

end
