function soleq = equilibrate(varargin)
% Uses initial conditions defined in DF and runs to equilibrium
%
%% Input arguments
% VARARGIN{1,1} = PAR
% VARARGIN{1,2} = ELECTRONIC_ONLY
% VARARGIN{1,3} = STORE_SOLUTION_WITH_EQM_TRAPS
% ELECTRONIC_ONLY:
% 0 = runs full equilibrate protocol
% 1 = skips ion equilibration
%
%% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
%% Start code
if length(varargin) == 1
    par = varargin{1,1};
    electronic_only = 0;
    store_eqm_traps = 0;
elseif length(varargin) == 2
    par = varargin{1,1};
    electronic_only = varargin{1,2};
    store_eqm_traps = 0;
elseif length(varargin) == 3
    par = varargin{1,1};
    electronic_only = varargin{1,2};
    store_eqm_traps = varargin{1,3};
else
    par = pc;
    electronic_only = 0;
    store_eqm_traps = 0;
end

tic;    % Start stopwatch
%% Initial arguments
% Setting sol.u = 0 enables a parameters structure to be read into
% DF but indicates that the initial conditions should be the
% analytical solutions
sol.u = 0;

% Store the original parameter set
par_origin = par;
% Start with zero SRH recombination
par.SRHset = 0;
% Radiative rec could initially be set to zero in addition if required
par.radset = 1;
% Start with no ionic carriers
par.eqm = 1;
% Switch off volumetric surface recombination check
par.vsr_check = 0;

%% General initial parameters
% Set applied bias to zero
par.V_fun_type = 'constant';
par.V_fun_arg(1) = 0;

% Set light intensities to zero
par.int1 = 0;
par.int2 = 0;
par.g1_fun_type = 'constant';
par.g2_fun_type = 'constant';

% Time mesh
par.tmesh_type = 2;
par.tpoints = 200;

% Series resistance
par.Rs = 0;

%% Switch off mobilities
par.mobset = 0;
par.mobseti = 0;
par.kineticset = 0;

%% Initial solution with zero mobility
disp('Initial solution, zero mobility')
sol = df(sol, par);
disp('Complete')

% Switch on mobilities
par.mobset = 1;
par.radset = 1;
par.SRHset = 1;

% Characteristic diffusion time
t_diff = (par.dcum0(end)^2)/(2*par.kB*par.T*min(min(par.mu_n), min(par.mu_p)));
par.tmax = 100*t_diff;
par.t0 = par.tmax/1e6;

%% Solution with mobility switched on
disp('Solution with mobility switched on')
sol = df(sol, par);

all_stable = verifyStabilization(sol.u, sol.t, 0.7);

% loop to check electrons have reached stable config
num = 0;
while any(all_stable) == 0 
    disp(['increasing equilibration time, tmax = ', num2str(par.tmax*10)]);
    par.tmax = 10*par.tmax;
    par.t0 = par.tmax/1e6;
    try
        sol = df(sol, par);
        all_stable = verifyStabilization(sol.u, sol.t, 0.7);
    catch
        warning('Stabilisation failed')
        num = num + 1;
        if num <= 2
            continue
        else
            break
        end
    end    
end
if store_eqm_traps
    soleq.elstatic = sol;
end
disp("Stabilisation verified");

Nt_guess = get_Nt_guess(sol);
sol.u(end,:,4) = Nt_guess;
if par_origin.taun >= par_origin.taup
    sol.par.taun(:) = 1;
    sol.par.taup(:) = 1*par_origin.taup./par_origin.taun;
else
    sol.par.taup(:) = 1;
    sol.par.taun(:) = 1*par_origin.taun./par_origin.taup;
end
sol.par.kineticset = 1;
sol.par.tmax = 100*t_diff;
sol.par.t0 = sol.par.tmax/1e3;
sol.par = refresh_device(sol.par);
sol.par.RelTol = 1e-6;
sol.par.AbsTol = 1e-9;
disp('Solution with mobility switched on and kinetic traps')
sol = df(sol);

all_stable = 0;
num = 0;
while any(all_stable) == 0 
    if num ~= 0
        disp(['increasing equilibration time, tmax = ', num2str(sol.par.tmax*10)]);
    end
    sol.par.tmax = 10*sol.par.tmax;
    sol.par.t0 = sol.par.tmax/1e6;
    try
        sol = df(sol);
        all_stable = verifyStabilization(sol.u(:,:,:), sol.t, 0.7);
        rec = dfana.calcr(sol, "sub") ;        
        rec_n = rec.srh_n(end,:);
        rec_p = rec.srh_p(end,:);
        [max_diff, arg_max_diff] = max(abs(rec_n-rec_p));
        max_rec = max(abs([rec_n(arg_max_diff),rec_p(arg_max_diff)]));
        if max_rec == 0
            max_rec = 1e-20;
        end
        if max_diff/max_rec > 1e-3 && max_rec > 1e2
            all_stable = 0.*all_stable;
        end
        num = num + 1;
    catch
        warning('Stabilisation failed')
        num = num + 1;
        if num <= 2
            continue
        else
            break
        end
    end    
end

disp("Stabilisation verified");

sol.par.taun(:) = par_origin.taun(:);
sol.par.taup(:) = par_origin.taup(:);
sol.par = refresh_device(sol.par);
soleq.el = sol;
%Make sure kinetic traps turned on
soleq.el.par.kineticset = 1;
soleq.el.par.N_ionic_species = 0;
soleq.el.par = refresh_device(soleq.el.par);

disp('Electronic carrier equilibration complete')

if electronic_only == 0 && par_origin.N_ionic_species > 0
    %% Equilibrium solutions with ion mobility switched on
    par.eqm = 0;

    % Create temporary solution for appending initial conditions to
    sol = soleq.el;
    sol.par.N_ionic_species = par_origin.N_ionic_species;

    % Start without SRH or series resistance
    %par.SRHset = 0;
    par.Rs = 0;

    % disp('Closed circuit equilibrium with ions')

    % Take ratio of electron and ion mobilities in the active layer
    rat_cation = par.mu_n(par.active_layer)./par.mu_ion(par.active_layer,:);
    
    for i = 1:length(rat_cation)
        if isnan(rat_cation(i)) || isinf(rat_cation(i))
            rat_cation(i) = 0;
        end
    end

    par.mobset = 1;
    par.mobseti = 1;           % Ions are accelerated to reach equilibrium
    par.K_ion = rat_cation';
    par.tmax = 1e4*t_diff;
    par.t0 = par.tmax/1e3;
    par.kineticset = 0;
    sol.par = par;

    disp('Closed circuit equilibrium with ions')

    sol = df(sol, par);
    all_stable = verifyStabilization(sol.u(:,:,:), sol.t, 0.7);

    % loop to check ions have reached stable config- if not accelerate ions by
    % order of mag
    num = 0;
    while any(all_stable) == 0 
        disp(['increasing equilibration time, tmax = ', num2str(par.tmax*10)]);
        par.tmax = par.tmax*10;
        par.t0 = par.tmax/1e6;
        try
            sol = df(sol, par);
            all_stable = verifyStabilization(sol.u(:,:,:), sol.t, 0.7);
        catch
            warning('Stabilisation failed')
            num = num + 1;
            if num <= 2
                continue
            else
                break
            end
        end
    end
    
    if store_eqm_traps
        soleq.ionstatic = sol;
        soleq.ionstatic.par.mobseti = 1;
        soleq.ionstatic.par.K_ion = ones(par.N_ionic_species,1);
    end
    disp("Stabilisation verified");
    
    Nt_guess = get_Nt_guess(sol);
    sol.u(end,:,4) = Nt_guess;
    sol.par.kineticset = 1;
    if par_origin.taun >= par_origin.taup
        sol.par.taun(:) = 1;
        sol.par.taup(:) = 1*par_origin.taup./par_origin.taun;
    else
        sol.par.taup(:) = 1;
        sol.par.taun(:) = 1*par_origin.taun./par_origin.taup;
    end
    sol.par.tmax = 1e4*t_diff;
    sol.par.t0 = sol.par.tmax/1e3;
    sol.par.K_ion = ones(par.N_ionic_species,1);
    sol.par = refresh_device(sol.par);

    disp('Closed circuit equilibrium with ions and kinetic traps')

    sol = df(sol);
    all_stable = 0;
    num = 0;
    while any(all_stable) == 0 
        disp(['increasing equilibration time, tmax = ', num2str(sol.par.tmax*10)]);
        sol.par.tmax = 10*sol.par.tmax;
        sol.par.t0 = sol.par.tmax/1e6;
        try
            sol = df(sol);
            all_stable = verifyStabilization(sol.u(:,:,:), sol.t, 0.7);
            rec = dfana.calcr(sol, "sub") ;
            rec_n = rec.srh_n(end,:);
            rec_p = rec.srh_p(end,:);
            [max_diff, arg_max_diff] = max(abs(rec_n-rec_p));
            max_rec = max(abs([rec_n(arg_max_diff),rec_p(arg_max_diff)])) ;
            if max_rec == 0
                max_rec = 1e-20;
            end
            if max_diff/max_rec > 1e-4 && max_rec > 1e2
                all_stable = 0.*all_stable;
            end
        catch
            warning('Stabilisation failed')
            num = num + 1;
            if num <= 2
                continue
            else
                break
            end
        end
    end

    disp("Stabilisation verified");
    
    sol.par.taun(:) = par_origin.taun(:);
    sol.par.taup(:) = par_origin.taup(:);
    sol.par = refresh_device(sol.par);
    % write solution
    soleq.ion = sol;
    soleq.ion.par.mobseti = 1;
    soleq.ion.par.kineticset = 1;
    soleq.ion.par.K_ion = ones(par.N_ionic_species,1);

    disp('Ionic carrier equilibration complete')
end

disp('EQUILIBRATION COMPLETE')
toc

end
