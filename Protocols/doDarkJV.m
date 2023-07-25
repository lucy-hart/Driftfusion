function J_dark = doDarkJV(sol_ini, voltage_ar, dwelltime_ar)
%Dark JV protocol
%Equilbirates the device in the dark for each voltage in voltage_ar for
%time in dwelltime_ar and then returns the current value at the final
%sampling point

% sol_ini   	= an initial solution
% voltage_ar    = voltage points to sample
% dwelltime_ar  = time to stabilise at each applied voltage

%%
%Code to deal with case where dwelltime_ar is just an integer i.e. all
%values have same dwell time
if length(dwelltime_ar) == 1
    dwelltime_ar = dwelltime_ar*ones(length(voltage_ar));
elseif length(dwelltime_ar) ~= length(voltage_ar)
    error('Dwell time must be an integer or an array of the same length as voltage_ar.')
end
%% General initial parameters
J_dark.sol = cell(1, length(voltage_ar));
par = sol_ini.par;

par.g1_fun_type = 'constant';
par.g1_fun_arg(1) = 0;        % For future proof
par.int1 = 0;

%%
%Find solution at each voltage in voltage_ar
for i = 1:length(voltage_ar)  
    disp(['Voltage = ' num2str(voltage_ar(i)) ' V'])
    %First sweep from 0 V to applied bias SLOWLY 
    par.tmesh_type = 1;
    par.t0 = 0;
    par.tmax = 100;

    par.V_fun_type = 'sweep';
    par.V_fun_arg(1) = 0;
    par.V_fun_arg(2) = voltage_ar(i);
    par.V_fun_arg(3) = par.tmax;

    sol = df(sol_ini, par);

    %Now hold at constant bias for 1 second
    par.tpoints = 500;
    par.tmax = dwelltime_ar(i);   
    
    par.V_fun_type = 'constant';
    par.V_fun_arg(1) = voltage_ar(i);
  
    J_dark.sol{i} = df(sol, par);
end
%%
% Get out the current values at each voltage
%Sample current at LH boundary (i.e. xpos = 0)
J_dark.Jvalue = zeros(1, length(voltage_ar));

figure('Name', 'J vs. t')
hold on 

for i = 1:length(voltage_ar)
    [~,t,~,~,~,~,~,~,~,~] = dfana.splitsol(J_dark.sol{i});
    [J, ~, ~] = dfana.calcJ(J_dark.sol{i});
    J_dark.Jvalue(i) = J.tot(end, 1);
    plot(t, J.tot(:,1), 'DisplayName', [num2str(voltage_ar(i)) ' V'])
end

hold off
xlabel('Time (s)')
ylabel('Current Density (mA cm^{-2})')
legend('Location', 'bestoutside')

        