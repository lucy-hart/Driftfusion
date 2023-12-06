function Vinvert = findVinvert_v2(JVsol)
%Finds the potentials at which the ion distribution invert on both sides of
%the device during a JV sweep
%For a symmetric device, these should be the same and so can use
%findVflation
%Not true for the asymmetric case and so need to find each seperately - in
%general device will be dped by injection from the contact layer with the
%better energetic alignment. Assume ETL for concreteness. Means cations
%accumulate in bulk of deivce, as well as at interface with the HTL and so
%the inversion votage at the HTL will be lower than the inversion voltage
%at the ETL
%Do this by finding where the sign of the point next to the interfaces
%changes sign (where sign refers to the sign of N(interface) - N_eqm
%%
epp = 8.85e-12;
kB = 1.38e-23;
e = 1.6e-19;

par = JVsol.par;
Vapp = dfana.calcVapp(JVsol);
maxV_point = ceil(length(Vapp)/2);

%Shift x coordiantes so 0 is the start of the perovskite layer
xstart = sum(par.layer_points(1:2)) + 1;
xstop = sum(par.layer_points(1:3));

%This is a num_voltage_values x num_active_layer_points matrix
ion_conc = JVsol.u(:,xstart:xstop+1,4) - par.Ncat(3);

ion_charge_HTL = ion_conc(:,1);
ion_charge_ETL = ion_conc(:,end);

if sign(ion_charge_HTL(1)) ~= sign(ion_charge_HTL(maxV_point))
    [~,argmin] = min(abs(ion_charge_HTL(1:maxV_point)));
    Vinvert.HTL = Vapp(argmin);
else
    Vinvert.HTL = 0;
end

if sign(ion_charge_ETL(1)) ~= sign(ion_charge_ETL(maxV_point))
    [~,argmin] = min(abs(ion_charge_ETL(1:maxV_point)));
    Vinvert.ETL = Vapp(argmin);
else
    Vinvert.ETL = 0;
end
