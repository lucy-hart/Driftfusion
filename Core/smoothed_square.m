function y = smoothed_square(coeff, t)
%Generates a smoothed out square wave to try and help solver converge

%% Start code
A_start = coeff(1);
A_pulse = coeff(2);
period = coeff(3);
duty_cycle = coeff(4);

t_pulse = (duty_cycle*period)/100;
tau = 1e-4;

t0 = 0.05*period;
t1 = t0 + t_pulse/2;

%lt = less than function
%returns an array with elements set to 1 when A<B; otherwise, the element
%is 0
%b = mod(a,m) returns the remainder after division of a by m (a modulo m)
%ge = greater than or equal to function
y = lt(t, t1).*(A_start + (A_pulse - A_start).*(1./(1+exp(-(t-t0)./tau)))) +...
    ge(t, t1).*(A_start + (A_pulse - A_start).*(1./(1+exp((t-t0-t_pulse)./tau)))); 

