params_filepath = './PEM_workshop_Input_files/1_layer_device_suns_Voc.csv';

%% Load in parameters
par = pc(params_filepath);
par.taun = 5e-8;
par.taup = 5e-8;
par = refresh_device(par);

%% Obtain equilibrium solution
sol_eq = equilibrate(par, 1);

Vmax = 1.3;                     % Maximum voltage for cyclic voltammogram
Vmin = -0.5;                    % Minimum voltage for cyclic voltammogram
scan_rate = 1e-3;               % Current-voltage scan rate [Vs-1] 

illuminations = logspace(-1,log10(5),15);
Vocs = zeros(1,length(illuminations));

for i = 1:length(illuminations)
    %% Call function to obtain equilibrium and cyclic voltammogram solutions
    sol_CV = doCV(sol_eq.el, illuminations(i), Vmin, Vmax, Vmin, scan_rate, 0.5, 181);
    Vapp = dfana.calcVapp(sol_CV);  % Applied voltage
    Jstruct = dfana.calcJ(sol_CV);
    J = Jstruct.tot(:,1)'; 
    Vocs(i) = interp1(J, Vapp, 0);
end
%% 
plot_suns_Voc(illuminations, Vocs, par)