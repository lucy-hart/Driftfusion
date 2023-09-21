function [V_Emin, Emin] = findEminbulk(JVsol)
%Finds the potential at which average electric field in the bulk is
%minimised during a JV sweep 

%Do this by calling dfana.calcF to get th E field and then fidning the
%voltage where the absolute value of this integraed across the ctive layer
%is minimised 

%%
par = JVsol.par;
Vapp = dfana.calcVapp(JVsol);

%Find indicaes of the x coordinates where the active layer starts and stops
xstart = sum(par.layer_points(1:2)) + 1;
xstop = sum(par.layer_points(1:3));
%This is a 1 x num_active_layer_points matrix
x = par.xx(xstart:xstop); 
%Calculate F
E = dfana.calcF(JVsol, "whole");
E_bulkonly = E(:,xstart:xstop);

integrated_abs_E = trapz(x, abs(E_bulkonly), 2);
plot(Vapp, integrated_abs_E)

[integrated_Emin, Eminarg] = min(integrated_abs_E);
%Divide integrated E by active layer width to get th expectaton value 
Emin = integrated_Emin/(x(end)-x(1));
V_Emin = Vapp(Eminarg);
end