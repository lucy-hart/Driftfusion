function r_srh_VBTS = calc_rsrh_VBTS(E_min, E_max, Evb, Ecb, E_UVB, Nt_VB, Cn_VB, Cp_VB, Nc, Nv, n, p, ni, kB, T)
% Calculates the conduction band SRH recombination rate between limits
% E_MIN and E_MAX
r_srh_VBTS_fun = @(E) (Nt_VB.*exp((Evb - E)./E_UVB)).*(Cn_VB.*Cp_VB.*(n*p - ni^2)./...
    (Cp_VB.*(p + Nv.*exp((Evb - E)./(kB*T))) + Cn_VB.*(n + Nc.*exp((E - Ecb)./(kB*T)))));
r_srh_VBTS = integral(r_srh_VBTS_fun, E_min, E_max);

end
