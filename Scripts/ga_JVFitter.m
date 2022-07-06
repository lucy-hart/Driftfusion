%Runs a ga on one of the devices used in Weidong_ETL. Aims to minimise the
%difference between the simulated and experimental JV parameters, as well as
% the SC and OC QFLS.

%% Read in data files 
par_kloc6 = pc('Input_files/PTAA_MAPI_Kloc6_v4.csv');
par_pcbm = pc('Input_files/PTAA_MAPI_PCBM_v4.csv');
par_icba = pc('Input_files/PTAA_MAPI_ICBA_v4.csv');
par_iph = pc('Input_files/PTAA_MAPI_IPH_v4.csv');

devices = {par_kloc6, par_pcbm, par_icba, par_iph};
which_device = 2;

%% Do the ga 

%Define bounds for varaibles
%Varying ETL LUMO, mobility and srv

UB = [1e-5 1e-3 1e-4 1e-3; 4.4 4.0 3.8 4.0; 1e4 1e3 1e3 1e3];
LB = [1e-7 1e-5 1e-6 1e-5; 4.2 3.9 3.7 3.8; 10 10 10 10];

% Using a function handle/ anaonymous function here. 
% Basically, @f(x)fun(x,y,z) tells MATLAB that, of the variables x, y 
% and z, x is the one which is being varied.
fitness = @(params)find_residuals(params, devices, which_device);

A = [];
b = [];
Aeq = [];
beq = [];
lb=LB(:,which_device)';
ub=UB(:,which_device)';

best_fit_params = ga(fitness, 3, A, b, Aeq, beq, lb, ub);

%%
function resid = find_residuals(params, devices, which_device)
%Calculates the residuals as a scalar value to pass to the ga

    %unpack variables from arguments passed by ga
    mob = params(1);
    LUMO = params(2);
    v_surf = params(3);
    
    Experimental = [19.8 23.0 15.6 21.2; 0.99 1.04 1.06 1.06; 0.46 0.79 0.30 0.78; 1.06 1.05 1.07 1.10; 1.07 1.10 1.12 1.13];
    Calculated = zeros(5,1);

    %find the residuals
    stats = test_JV(devices{which_device}, 1, mob, LUMO, v_surf);
    if isa(stats,'struct') == 1
        Calculated(1,1) = abs(stats.Jsc_f);
        Calculated(2,1) = stats.Voc_f;
        Calculated(3,1) = stats.FF_f;
        Calculated(4,1) = stats.QFLS_SC;
        Calculated(5,1) = stats.QFLS_OC;
    elseif isa(stats,'struct') == 0
        Calculated(1,1) = 0;
        Calculated(2,1) = 0;
        Calculated(3,1) = 0;
        Calculated(4,1) = 0;
        Calculated(5,1) = 0;
    end 

    Residual_matrix = ((Experimental(:,which_device)-Calculated)./Experimental(:,which_device)).^2; 
    resid = sum(Residual_matrix);

end

