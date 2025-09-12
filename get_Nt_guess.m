function Nt_guess = get_Nt_guess(sol)
    dev = sol.par.dev;

    krec_n = 1./(dev.taun.*dev.Ntrap_n);
    krec_p = 1./(dev.taup.*dev.Ntrap_n);
    kout_p = dev.pt.*krec_p;

    n = sol.u(end,:,2);
    p = sol.u(end,:,3);

    Nt_guess = dev.Ntrap_n.*((krec_n.*n+kout_p)./(krec_n.*(p+dev.pt)+krec_p.*(n+dev.nt)));

end