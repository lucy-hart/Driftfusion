function ptrap_VBTS = calc_ptrap_VBTS(E_min, E_max, Evb, Ecb, E_UVB, Nt_VB, Cn_VB, Cp_VB, Nc, Nv, n, p, ni, kB, T)
% Calculates the valence band tail state trapped hole density between limits
% E_MIN and E_MAX
ptrap_VBTS_fun = @(E) (Nt_VB.*exp((Evb - E)./E_UVB)).*((Cp_VB*p + Cn_VB.*Nc.*exp((E - Ecb)./(kB*T)))./...
    (Cn_VB.*(n + Nc.*exp((E - Ecb)./(kB*T))) + Cp_VB.*(p + Nv.*exp((Evb - E)./(kB*T)))));

ptrap_VBTS = integral(ptrap_VBTS_fun, E_min, E_max);

end
