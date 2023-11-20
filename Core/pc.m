classdef pc
% PC (Parameters Class) defines all the required properties for your
% device. PC.BUILDDEV builds a structure PO.DEV (where PO is a Parameters Object)
% that defines the properties of the device at every spatial mesh point, including
% interfaces. Whenever PROPERTIES are overwritten in a protocol, the device should
% be rebuilt manually using PC.BUILDDEV. The spatial mesh is a linear piece-wise mesh
% and is built by the MESHGEN_X function. Details of how to define the mesh
% are given below in the SPATIAL MESH SUBSECTION.
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
    properties (Constant)
        %% Physical constants
        kB = 8.617330350e-5;     % Boltzmann constant [eV K^-1]
        epp0 = 552434;           % Epsilon_0 [e^2 eV^-1 cm^-1] - Checked (02-11-15)
        q = 1;                   % Charge of the species in units of e.
        e = 1.60217662e-19;      % Elementary charge in Coulombs.
    end

    properties
        % Temperature [K]
        T = 300;

        %% Spatial mesh
        % Device Dimensions [cm]
        % The spatial mesh is a linear piece-wise mesh and is built by the
        % MESHGEN_X function using 2 arrays DCELL and PCELL,
        % which define the thickness and number of points of each layer
        % respectively.
        d = 400e-7;         % Layer and subsection thickness array
        layer_points = 400;            % Points array

        %% Layer description
        % Define the layer type for each of the layers in the device. The
        % options are:
        % LAYER = standard layer
        % ACTIVE = standard layer but the properties of this layer are
        % flagged such that they can easily be accessed
        % JUNCTION = a region with graded properties between two materials
        % (either LAYER or ACTIVE type)
        % with different properties
        layer_type = {'active'}
        % STACK is used for reading the optical properties library.
        % The names here do not influence the electrical properties of the
        % device. See INDEX OF REFRACTION LIBRARY for choices- names must be entered
        % exactly as given in the column headings with the '_n', '_k' omitted
        material = {'MAPICl'}
        layer_colour = [1,1,1;1,1,1;1,1,1;1,1,1;1,1,1;1,1,1];
        % Define spatial cordinate system- typically this will be kept at
        % 0 for most applications
        % m=0 cartesian
        % m=1 cylindrical polar coordinates
        % m=2 spherical polar coordinates
        m = 0;

        %% Spatial mesh
        % xmesh_type specification - see MESHGEN_X.
        xmesh_type = 'erf-linear';
        xmesh_coeff = [0.7];        % Coefficient array for defining point spacing
        
        %% Time mesh
        % The time mesh is dynamically generated by ODE15s- the mesh
        % defined by MESHGEN_T only defines the values of the points that
        % are read out and so does not influence convergence. Defining an
        % unecessarily high number of points however can be expensive owing
        % to interpolation of the solution.
        tmesh_type = 'log10';             % Mesh type- for use with meshgen_t
        t0 = 1e-16;                 % Initial log mesh time value
        tmax = 1e-12;               % Max time value
        tpoints = 100;              % Number of time points

        %% GENERAL CONTROL PARAMETERS
        mobset = 1;                         % Switch on/off electron hole mobility- MUST BE SET TO ZERO FOR INITIAL SOLUTION
        mobseti = 1;                        % Switch on/off ionic carrier mobility- MUST BE SET TO ZERO FOR INITIAL SOLUTION
        SRHset = 1;                         % Switch on/off SRH recombination - recommend setting to zero for initial solution
        radset = 1;                         % Switch on/off band-to-band recombination
        N_max_variables = 5;                % Total number of allowable variables in this version of Driftfusion
        prob_distro_function = 'Blakemore'; % 'Fermi' = Fermi-Dirac, 'Blakemore' = Blakemore aproximation, 'Boltz' = Boltzmann statistics
        gamma_Blakemore = 0.27;                       % Blakemore coefficient    
        Fermi_limit = 0.2;                  % Max allowable limit for Fermi levels beyond the bands [eV]
        Fermi_Dn_points = 400;              % No. of points in the Fermi-Dirac look-up table
        intgradfun = 'linear'               % Interface gradient function 'linear' = linear, 'erf' = 'error function'

        %% Generation
        % optical_model = Optical Model
        % 0 = Uniform Generation
        % 1 = Beer Lambert
        optical_model = 'Beer-Lambert';
        int1 = 0;               % Light intensity source 1 (multiples of g0 or 1 sun for Beer-Lambert)
        int2 = 0;               % Light intensity source 2 (multiples of g0 or 1 sun for Beer-Lambert)
        g0 = [2.6409e+21];      % Uniform generation rate [cm-3s-1]
        light_source1 = 'AM15';
        light_source2 = 'laser';
        laser_lambda1 = 0;
        laser_lambda2 = 638;
        g1_fun_type = 'constant'
        g2_fun_type = 'constant'
        g1_fun_arg = 0;
        g2_fun_arg = 0;
        side = 'left';                           % illumination side 1 = left, 2 = right
        % default: Approximate Uniform generation rate @ 1 Sun for 510 nm active layer thickness

        %% Pulse settings
        pulsepow = 10;          % Pulse power [mW cm-2] OM2 (Beer-Lambert and Transfer Matrix only)

        %%%%%%%%%%% LAYER MATERIAL PROPERTIES %%%%%%%%%%%%%%%%%%%%
        % Numerical values should be given as a row vector with the number of
        % entries equal to the number of layers specified in STACK

        %% Energy levels [eV]
        Phi_EA = [0];           % Electron affinity
        Phi_IP = [-1];           % Ionisation potential

        %% Equilibrium Fermi energies [eV]
        % These define the doping density in each layer- see NA and ND calculations in methods
        EF0 = [-0.5];

        %% SRH trap energies [eV]
        % These must exist within the energy gap of the appropriate layers
        % and define the variables PT and NT in the expression:
        % U = (np-ni^2)/(taun(p+pt) +taup(n+nt))
        Et =[-0.5];
        Et2 =[-0.5];
        ni_eff = 0;     % Effective intrinsic carrier density used for surface recombination equivalence
        
        %% Electrode Fermi energies [eV]
        % Fermi energies of the metal electrode. These define the built-in voltage, Vbi
        % and the boundary carrier concentrations n0_l, p0_l, n0_r, and
        % p0_r
        Phi_left = -0.6;
        Phi_right = -0.4;

        %% Effective Density Of States (eDOS) [cm-3]
        Nc = [1e19];
        Nv = [1e19];
        % PEDOT eDOS: https://aip.scitation.org/doi/10.1063/1.4824104
        % MAPI eDOS: F. Brivio, K. T. Butler, A. Walsh and M. van Schilfgaarde, Phys. Rev. B, 2014, 89, 155204.
        % PCBM eDOS:
        
        %% Mobile ions        
        N_ionic_species = 1;        
        Nani = [1e19];                  % Mobile ion defect density [cm-3] - A. Walsh et al. Angewandte Chemie, 2015, 127, 1811.
        Ncat = [1e19];                  % Mobile ion defect density [cm-3] - A. Walsh et al. Angewandte Chemie, 2015, 127, 1811.
        z_c = 1;                        % Integer charge state for cations
        z_a = -1;                       % Integer charge state for anions
        % Limits the density of ions - Approximate density of iodide sites [cm-3]
        a_max = [1.21e22];                 % P. Calado thesis
        c_max = [1.21e22];
        
        K_a = 1;                    % Coefficients to easily accelerate ions
        K_c = 1;                   % Coefficients to easily accelerate ions
        
        %% Mobilities   [cm2V-1s-1]
        mu_n = [1];         % electron mobility
        mu_p = [1];         % hole mobility
        mu_c = [1e-10];
        mu_a = [1e-12]; 
        % PTPD h+ mobility: https://pubs.rsc.org/en/content/articlehtml/2014/ra/c4ra05564k
        % PEDOT mu_n = 0.01 cm2V-1s-1 https://aip.scitation.org/doi/10.1063/1.4824104
        % TiO2 mu_n = 0.09 cm2V-1s-1 Bak2008
        % Spiro mu_p = 0.02 cm2V-1s-1 Hawash2018
        %% Relative dielectric constants
        epp = [10];
        epp_factor = 1e6;    % a factor required to a prevent singular matrix- still under investigation
        %% Recombination
        % Radiative recombination, r_rad = k(np - ni^2)
        % [cm3 s-1] Radiative Recombination coefficient
        B = [3.6e-12];

        %% SRH time constants for each layer [s]
        taun = [1e6];           % [s] SRH time constant for electrons
        taup = [1e6];           % [s] SRH time constant for holes  
        
        %% Surface recombination and extraction coefficients [cm s-1]
        % Descriptions given in the comments considering that holes are
        % extracted at left boundary, electrons at right boundary
        sn_l = 1e7;     % electron surface recombination velocity left boundary
        sn_r = 1e7;     % electron extraction velocity right boundary
        sp_l = 1e7;     % hole extraction left boundary
        sp_r = 1e7;     % hole surface recombination velocity right boundary
        
        %% Volumetric surface recombination
        vsr_mode = 1;               % Either 1 for volumetric surface recombination approximation or 0 for off
        vsr_check = 1;              % Perform check for self-consitency at the end of DF
        sn = [0];                   % Electron interfacial surface recombination velocity [cm s-1]
        sp = [0];                   % Hole interfacial surface recombination velocities [cm s-1]
        sn2 = [0];                  % Electron interfacial surface recombination velocity [cm s-1]
        sp2 = [0];                  % Hole interfacial surface recombination velocities [cm s-1]
        frac_vsr_zone = 0.1;        % recombination zone thickness [fraction of interface thickness]
        vsr_zone_loc = {'auto'};    % recombination zone location either: 'L', 'C', 'R', or 'auto'. IMPORT_PROPERTIES deals with the choice of value.
        AbsTol_vsr = 1e10;          % The integrated interfacial recombination flux above which a warning can be flagged [cm-2 s-1]
        RelTol_vsr = 0.05;          % Fractional error between abrupt and volumetric surface recombination models above which a warning is flagged
        
        %% Series resistance
        Rs = 0;
        Rs_initial = 0;         % Switch to allow linear ramp of Rs on first application

        %% Defect recombination rate coefficient
        % Currently not used
        k_defect_p = 0;
        k_defect_n = 0; 

        %% Dynamically created variables
        genspace = [];
        x = [];
        xx = [];
        x_sub = [];
        dev = [];
        dev_sub = [];
        t = [];
        xpoints = [];
        gx1 = [];       % Light source 1
        gx2 = [];       % Light source 2

        %% Voltage function parameters
        V_fun_type = 'constant';
        V_fun_arg = 0;
        
        % Define the default relative tolerance for the pdepe solver
        % 1e-3 is the default, can be decreased if more precision is needed
        % Solver options
        MaxStepFactor = 1;      % Multiplier for easy access to maximum time step
        RelTol = 1e-3;
        AbsTol = 1e-6;
        
        %% Impedance parameters
        J_E_func = [];
        J_E_func_tilted = [];
        E2_func = [];
    end

    %%  Properties whose values depend on other properties (see 'get' methods).
    properties (Dependent)
        active_layer
        dcell
        parr
        d_active
        dcum
        dcum0           % includes first entry as zero
        d_midactive
        dEAdx
        dIPdx
        dNcdx
        dNvdx
        gamma
        int_switch
        Dn
        Eg
        Efi
        NA
        ND
        Vbi
        n0
        n0_l
        n0_r
        ni
        nt              % Density of CB electrons when Fermi level at trap state energy
        nt2
        nt_inter
        p0
        pcum
        pcum0           % Includes first entry as zero
        p0_l
        p0_r
        pt              % Density of VB holes when Fermi level at trap state energy
        pt2
        pt_inter
        wn
        wp
        wscr            % Space charge region width
        x0              % Initial spatial mesh value
    end

    methods
        function par = pc(varargin)
            % Parameters constructor function- runs numerous checks that
            % the input properties are consistent with the model
            if length(varargin) == 1
                % Use argument as filepath and overwrite properties using
                % PC.IMPORTPROPERTIES
                filepath = varargin;
                par = import_properties(par, filepath);
            elseif length(varargin) > 1
                filepath = varargin{1, 1};
                par = import_properties(par, filepath);
                warning('pc should have 0 or 1 input arguments- only the first argument will be used for the filepath')
            end

            % Warn if xmesh_type is not correct
            if ~ any(strcmp(par.xmesh_type, {'linear', 'erf-linear'}))
                error('PAR.xmesh_type should either be ''linear'' or ''erf-linear''. MESHGEN_X cannot generate a mesh if this is not the case.')
            end

            % Warn if doping density exceeds eDOS
            for i = 1:length(par.ND)
                if par.ND(i) > par.Nc(i) || par.NA(i) > par.Nc(i)
                    msg = 'Doping density must be less than eDOS. For consistent values ensure electrode workfunctions are within the band gap and check expressions for doping density in Dependent variables.';
                    error(msg);
                end
            end

            % Warn if trap energies are outside of band gap energies
            for i = 1:length(par.Et)
                if par.Et(i) >= par.Phi_EA(i) || par.Et(i) <= par.Phi_IP(i)
                    msg = 'Trap energies must exist within layer band gap.';
                    error(msg);
                end
            end

            % Warn if a_max is set to zero in any layers - leads to
            % infinite diffusion rate
            for i = 1:length(par.a_max)
                if par.a_max(i) <= 0
                    msg = 'Maximum cation density (a_max) cannot have zero or negative entries- choose a low value rather than zero e.g. 1';
                    error(msg);
                end
            end

            % Warn if c_max is set to zero in any layers - leads to
            % infinite diffusion rate
            for i = 1:length(par.c_max)
                if par.c_max(i) <= 0
                    msg = 'Maximum cation density (c_max) cannot have zero or negative entries- choose a low value rather than zero e.g. 1';
                    error(msg);
                end
            end
            
            % Warn if electrode workfunctions are outside of boundary layer
            % bandgap
            if par.Phi_left < par.Phi_IP(1) || par.Phi_left > par.Phi_EA(1)
                msg = 'Left-hand workfunction (Phi_left) out of range: value must exist within left-hand layer band gap';
                error(msg)
            end

            if par.Phi_right < par.Phi_IP(end) || par.Phi_right > par.Phi_EA(end)
                msg = 'Right-hand workfunction (Phi_right) out of range: value must exist within right-hand layer band gap';
                error(msg)
            end

            % Warn if property array do not have the correct number of
            % layers. The layer thickness array is used to define the
            % number of layers
            if length(par.parr) ~= length(par.d)
                msg = 'Points array (parr) does not have the correct number of elements. Property arrays must have the same number of elements as the thickness array (d), except SRH properties for interfaces which should have length(d)-1 elements.';
                error(msg);
            elseif length(par.Phi_EA) ~= length(par.d)
                msg = 'Electron Affinity array (Phi_EA) does not have the correct number of elements. Property arrays must have the same number of elements as the thickness array (d), except SRH properties for interfaces which should have length(d)-1 elements.';
                error(msg);
            elseif length(par.Phi_IP) ~= length(par.d)
                msg = 'Ionisation Potential array (Phi_IP) does not have the correct number of elements. Property arrays must have the same number of elements as the thickness array (d), except SRH properties for interfaces which should have length(d)-1 elements.';
                error(msg);
            elseif length(par.mu_n) ~= length(par.d)
                msg = 'Electron mobility array (mu_n) does not have the correct number of elements. Property arrays must have the same number of elements as the thickness array (d), except SRH properties for interfaces which should have length(d)-1 elements.';
                error(msg);
            elseif length(par.mu_p) ~= length(par.d)
                msg = 'Hole mobility array (mu_n) does not have the correct number of elements. Property arrays must have the same number of elements as the thickness array (d), except SRH properties for interfaces which should have length(d)-1 elements.';
                error(msg);
            elseif length(par.mu_a) ~= length(par.d)
                msg = 'Ion mobility array (mu_p) does not have the correct number of elements. Property arrays must have the same number of elements as the thickness array (d), except SRH properties for interfaces which should have length(d)-1 elements.';
                error(msg);
            elseif length(par.NA) ~= length(par.d)
                msg = 'Acceptor density array (NA) does not have the correct number of elements. Property arrays must have the same number of elements as the thickness array (d), except SRH properties for interfaces which should have length(d)-1 elements.';
                error(msg);
            elseif length(par.ND) ~= length(par.d)
                msg = 'Donor density array (ND) does not have the correct number of elements. Property arrays must have the same number of elements as the thickness array (d), except SRH properties for interfaces which should have length(d)-1 elements.';
                error(msg);
            elseif length(par.Nc) ~= length(par.d)
                msg = 'Effective density of states array (Nc) does not have the correct number of elements. Property arrays must have the same number of elements as the thickness array (d), except SRH properties for interfaces which should have length(d)-1 elements.';
                error(msg);
            elseif length(par.Nv) ~= length(par.d)
                msg = 'Effective density of states array (Nv) does not have the correct number of elements. Property arrays must have the same number of elements as the thickness array (d), except SRH properties for interfaces which should have length(d)-1 elements.';
                error(msg);
            elseif length(par.Nani) ~= length(par.d)
                msg = 'Background ion density (Nani) does not have the correct number of elements. Property arrays must have the same number of elements as the thickness array (d), except SRH properties for interfaces which should have length(d)-1 elements.';
                error(msg);
            elseif length(par.a_max) ~= length(par.d)
                msg = 'Ion density of states array (a_max) does not have the correct number of elements. Property arrays must have the same number of elements as the thickness array (d), except SRH properties for interfaces which should have length(d)-1 elements.';
                error(msg);
            elseif length(par.epp) ~= length(par.d)
                msg = 'Relative dielectric constant array (epp) does not have the correct number of elements. Property arrays must have the same number of elements as the thickness array (d), except SRH properties for interfaces which should have length(d)-1 elements.';
                error(msg);
            elseif length(par.B) ~= length(par.d)
                msg = 'Radiative recombination coefficient array (B) does not have the correct number of elements. Property arrays must have the same number of elements as the thickness array (d), except SRH properties for interfaces which should have length(d)-1 elements.';
                error(msg);
            elseif length(par.EF0) ~= length(par.d)
                msg = 'Equilibrium Fermi level array (EF0) does not have the correct number of elements. Property arrays must have the same number of elements as the thickness array (d), except SRH properties for interfaces which should have length(d)-1 elements.';
                error(msg);
            elseif length(par.g0) ~= length(par.d)
                msg = 'Uniform generation array (g0) does not have the correct number of elements. Property arrays must have the same number of elements as the thickness array (d), except SRH properties for interfaces which should have length(d)-1 elements.';
                error(msg);
            elseif length(par.taun) ~= length(par.d)
                msg = 'Bulk SRH electron time constants array (taun_bulk) does not have the correct number of elements. Property arrays must have the same number of elements as the thickness array (d), except SRH properties for interfaces which should have length(d)-1 elements.';
                error(msg);
            elseif length(par.taup) ~= length(par.d)
                msg = 'Bulk SRH hole time constants array (taup_bulk) does not have the correct number of elements. Property arrays must have the same number of elements as the thickness array (d), except SRH properties for interfaces which should have length(d)-1 elements.';
                error(msg);
            elseif length(par.Et) ~= length(par.d)
                msg = 'Bulk SRH trap energy array (Et) does not have the correct number of elements. Property arrays must have the same number of elements as the thickness array (d), except SRH properties for interfaces which should have length(d)-1 elements.';
                error(msg);
            end

            %% Device and generation builder
            % Import variables and structure, xx, gx1, gx2, and dev must be
            % refreshed when to rebuild the device for example when
            % changing device thickness on the fly. These are not present
            % in the dependent variables as it is too costly to have them
            % continuously called.
            par = refresh_device(par);
        end
                
        function par = set.xmesh_type(par, value)
            if isa(value, 'double')
                % Backwards compat values
                switch value
                    case 4
                        par.xmesh_type = 'linear';
                    case 5
                        par.xmesh_type = 'erf-linear';
                end
            elseif isa(value, 'cell')
                if any(strcmp(par.xmesh_type, {'linear', 'erf-linear'}))
                    par.xmesh_type = value{1};
                end
            elseif isa(value, 'char')
                if any(strcmp(par.xmesh_type, {'linear', 'erf-linear'}))
                    par.xmesh_type = value;
                end
            else
                par.xmesh_type = 'erf-linear';
                warning('xmesh_type not recognised- defaulting to ''erf-linear'' spatial mesh');
            end
        end

        function par = set.tmesh_type(par, value)
            % Backwards compat values
            if isa(value, 'double')
                switch value
                    case 1
                        par.tmesh_type = 'linear';
                    case 2
                        par.tmesh_type = 'log10';
                end
            elseif isa(value, 'cell')
                if any(strcmp(value, {'linear', 'log10', 'log10-double'}))
                    par.tmesh_type = value{1};
                end
            elseif isa(value, 'char')
                if any(strcmp(value, {'linear', 'log10', 'log10-double'}))
                    par.tmesh_type = value;
                end
            else
                par.tmesh_type = 'linear';
                warning('tmesh_type not recognised- defaulting to ''linear'' mesh');
            end
        end
        
        function par = set.optical_model(par, value)
            % Backwards compat values
            if isa(value, 'double')
                switch value
                    case 0
                        par.optical_model = 'uniform';
                    case 1
                        par.optical_model = 'Beer-Lambert';
                end
            elseif isa(value, 'cell')
                if any(strcmp(value, {'uniform', 'Beer-Lambert'}))
                    par.optical_model = value{1};
                end
            elseif isa(value, 'char')
                if any(strcmp(value, {'uniform', 'Beer-Lambert'}))
                    par.optical_model = value;
                end
            else
                par.optical_model = 'Beer-Lambert';
                warning('optical_model not recognised- defaulting to ''Beer-Lambert''');
            end
        end
        
        function par = set.side(par, value)
            % Backwards compat values
            if isa(value, 'double')
                switch value
                    case 1
                        par.side = 'left';
                    case 2
                        par.side = 'right';
                end
            elseif isa(value, 'cell')
                if any(strcmp(value, {'left', 'right'}))
                    par.side = value{1};
                end
            elseif isa(value, 'char')
                if any(strcmp(value, {'left', 'right'}))
                    par.side = value;
                end
            else
                par.side = 'left';
                warning('illumination side not recognised- defaulting to ''left''');
            end
        end
        
        function par = set.taun(par, value)
            for i = 1:length(value)
                if isnan(value(i))
                    par.taun(i) = 1e100;
                else
                    par.taun(i) = value(i);
                end
            end
        end

        function par = set.taup(par, value)
            for i = 1:length(value)
                if isnan(value(i))
                    par.taup(i) = 1e100;
                else
                    par.taup(i) = value(i);
                end
            end
        end
        
        function par = set.ND(par, value)
            for i = 1:length(par.ND)
                if value(i) >= par.Nc(i)
                    error('Doping density must be less than eDOS. For consistent values ensure electrode workfunctions are within the band gap.')
                end
            end
        end

        function par = set.NA(par, value)
            for i = 1:length(par.ND)
                if value(i) >= par.Nv(i)
                    error('Doping density must be less than eDOS. For consistent values ensure electrode workfunctions are within the band gap.')
                end
            end
        end

        function value = get.gamma(par)
            switch par.prob_distro_function
                case 'Boltz'
                    value = 0;
                case 'Blakemore'
                    value = par.gamma_Blakemore;
            end
        end
        %% Get active layer indexes from layer_type
        function value = get.active_layer(par)
            value = find(strncmp('active', par.layer_type, 6));
            if length(value) == 0
                % If no flag is give assume active layer is middle
                value = round(length(par.layer_type)/2);
                warning('No designated ''active'' layer- assigning middle layer to be active')
            end
        end
        
        %% Layer thicknesses [cm]
        function value = get.dcell(par)
            % For backwards comptibility. layer_points and parr arre the now the
            % same thing
            value = par.d;
        end

        %% Layer points
        function value = get.parr(par)
            % For backwards comptibility. layer_points and parr arre the now the
            % same thing
            value = par.layer_points;
        end

        %% Active layer thickness
        function value = get.d_active(par)
            value = sum(par.dcell(par.active_layer(1):par.active_layer(end)));
        end
    
        function value = get.d_midactive(par)
           value = par.dcum0(par.active_layer(1)) + par.d_active/2;
        end
        
        %% Band gap energies    [eV]
        function value = get.Eg(par)
            value = par.Phi_EA - par.Phi_IP;
        end

        %% Built-in voltage Vbi based on difference in boundary workfunctions
        function value = get.Vbi(par)
            value = par.Phi_right - par.Phi_left;
        end

        %% Intrinsic Fermi Energies
        % Currently uses Boltzmann stats as approximation should always be
        function value = get.Efi(par)
            value = 0.5.*(par.Phi_EA+par.Phi_IP)+par.kB*par.T*log(par.Nc./par.Nv);
        end

        %% Donor densities
        function value = get.ND(par)
            value = distro_fun.nfun(par.Nc, par.Phi_EA, par.EF0, par);
        end

        %% Acceptor densities
        function value = get.NA(par)
            value = distro_fun.pfun(par.Nv, par.Phi_IP, par.EF0, par);
        end
        
        %% Intrinsic carrier densities (Boltzmann)
        function value = get.ni(par)
            value = ((par.Nc.*par.Nv).^0.5).*exp(-par.Eg./(2*par.kB*par.T));
        end

        %% Equilibrium electron densities
        function value = get.n0(par)
            value = distro_fun.nfun(par.Nc, par.Phi_EA, par.EF0, par);
        end

        %% Equilibrium hole densities
        function value = get.p0(par)
            value = distro_fun.pfun(par.Nv, par.Phi_IP, par.EF0, par);

        end

        %% Boundary electron and hole densities
        % Uses metal Fermi energies to calculate boundary densities
        % Electrons left boundary
        function value = get.n0_l(par)
            value = distro_fun.nfun(par.Nc(1), par.Phi_EA(1), par.Phi_left, par);
        end

        % Electrons right boundary
        function value = get.n0_r(par)
            value = distro_fun.nfun(par.Nc(end), par.Phi_EA(end), par.Phi_right, par);
        end

        % Holes left boundary
        function value = get.p0_l(par)
            value = distro_fun.pfun(par.Nv(1), par.Phi_IP(1), par.Phi_left, par);
        end

        % holes right boundary
        function value = get.p0_r(par)
            value = distro_fun.pfun(par.Nv(end), par.Phi_IP(end), par.Phi_right, par);
        end
       
        %% SRH trap energy coefficients
        function value = get.nt(par)
            value = distro_fun.nfun(par.Nc, par.Phi_EA, par.Et, par);
        end
        
        function value = get.pt(par)
            value = distro_fun.pfun(par.Nv, par.Phi_IP, par.Et, par);
        end

        function value = get.nt2(par)
            value = distro_fun.nfun(par.Nc, par.Phi_EA, par.Et2, par);
        end
        
        function value = get.pt2(par)
            value = distro_fun.pfun(par.Nv, par.Phi_IP, par.Et2, par);
        end
        
        %% Thickness and point arrays
        function value = get.dcum(par)
            value = cumsum(par.dcell);
        end

        function value = get.pcum(par)
            value = cumsum(par.layer_points);
        end

        function value = get.pcum0(par)
            value = [1, cumsum(par.layer_points)];
        end

        function value = get.dcum0(par)
            value = [0, cumsum(par.dcell)];
        end
        
        % interface switch for zeroing field in interfaces
        function value = get.int_switch(par)
            value = ones(1, length(par.material));
        end
        
    end

    methods (Static)

        function xx = xmeshini(par) % For backwards compatibility
            xx = meshgen_x(par);
        end     

    end
end
