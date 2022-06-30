function ntrap_CBTS = calc_ntrap_CBTS(E_min, E_max, Evb, Ecb, E_UCB, Nt_CB, Cn_CB, Cp_CB, Nc, Nv, n, p, ni, kB, T)
% Calculates the conduction band tail state trapped electron density between limits
% E_MIN and E_MAX

ntrap_CBTS_fun = @(E) (Nt_CB.*exp((E - Ecb)./E_UCB)).*((Cn_CB*n + Cp_CB*Nv.*exp((Evb - E)./(kB*T)))./...
    (Cn_CB.*(n + Nc.*exp((E - Ecb)./(kB*T))) + Cp_CB.*(p + Nv.*exp((Evb - E)./(kB*T)))));
ntrap_CBTS = integral(ntrap_CBTS_fun, E_min, E_max);

end
