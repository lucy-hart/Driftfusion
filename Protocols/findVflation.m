function Vflat = findVflation(JVsol)
%Finds the potential at which the ion distribution is most flat during a JV
%sweep. 
%Do this by calculating dipole moment of ionic charge in the active layer
%and interpolating to find where the voltage is at a minimum

%%
par = JVsol.par;
Vapp = dfana.calcVapp(JVsol);

%Shift x coordiantes so centre is at the centre of the active layer
%Do this as I want to find the dipole moment of the ionic charge
xstart = sum(par.layer_points(1:2)) + 1;
xstop = sum(par.layer_points(1:3));
%This is a 1 x num_active_layer_points matrix
x = par.xx(xstart:xstop) - (par.xx(1)+par.xx(end))/2;
%This is a num_voltage_values x num_active_layer_points matrix
ion_conc = JVsol.u(:,xstart:xstop,4) - par.Ncat(3);

ion_dipole_moment = trapz(x, x.*ion_conc, 2);

Vflat = interp1(ion_dipole_moment, Vapp, 0);
end