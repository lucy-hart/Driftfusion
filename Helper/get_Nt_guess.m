function Nt_guess = get_Nt_guess(varargin)
    if length(varargin) == 1
        sol = varargin{1,1};
        i = length(sol.t);
    elseif length(varargin) == 2
        sol = varargin{1,1};
        i = varargin{1,2};
    end
    dev = sol.par.dev;

    krec_n = 1./(dev.taun.*dev.Ntrap_n);
    krec_p = 1./(dev.taup.*dev.Ntrap_n);
    kout_p = dev.pt.*krec_p;

    n = sol.u(i,:,2);
    p = sol.u(i,:,3);

    Nt_guess = dev.Ntrap_n.*((krec_n.*n+kout_p)./(krec_n.*(p+dev.pt)+krec_p.*(n+dev.nt)));

end