%Code to calculate charge stored on device active layers
sol = eqm_QJV.el;

num_start = sum(sol.par.layer_points(1:2))+1;
num_stop = num_start + sol.par.layer_points(3)-1;
x = sol.x;
q = par.e*trapz(x(num_start:num_stop), sol.u(:,num_start:num_stop,2),2);

% num_start_ion = sol.par.layer_points(1)/2;
% num_stop_ion = num_start_ion + sol.par.layer_points(1)/2;
% qion = par.e*trapz(x(num_start_ion:num_stop_ion), sol.u(:,num_start_ion:num_stop_ion,4) - par.Ncat(2),2);
% qion2 = par.e*trapz(x(1:num_start_ion), sol.u(:,1:num_start_ion,4) - par.Ncat(2),2);





