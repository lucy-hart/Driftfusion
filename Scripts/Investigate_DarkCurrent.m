% parC60 = pc('Input_files/HTL_MAPI_C60_DavideValues.csv');
%The _dark file has a finer point sapcing as I wanted to see if this
%reduced the zero error
%It did not...
% parC60 = pc('Input_files/HTL_MAPI_C60_DavideValues_dark.csv');
% parPM6 = pc('Input_files/HTL_MAPI_PM6Y6_DavideValues.csv');

% eqm_C60 = equilibrate(parC60);
% eqm_PM6 = equilibrate(parPM6);

%%
%Test using pn junction
test_pn = 1;
if test_pn == 1
    par_pn_junction = pc('Input_files/pn_junction.csv');
    par_pn_junction.prob_distro_function = 'Boltz';
    par_pn_junction = refresh_device(par_pn_junction);
    eqm_pn = equilibrate(par_pn_junction);

    par = eqm_pn.el.par;
%Calculate saturation current using expression from https://en.wikipedia.org/wiki/Saturation_current
%First layer is the p-type and second layer is the n-type
%Various factors of 10 to the power of stuff are to convert everything to
%metres and then I've divided by 1e4 to get units of A cm-2
%par.kB returns kB in units of eV K-1, whixh is what you need to use the
%Einstien relation (see notes, 28/07/2)
    J_sat_calc = 1e2*par.e*par.ni(1)^2*((sqrt(par.kB*par.T*1e-4*par.mu_p(4)/par.taup(4))/par.ND(4)) + (sqrt(par.kB*par.T*1e-4*par.mu_n(1)/par.taun(1))/par.NA(1)));

    voltage_ar = linspace(-0.5, 0.2, 15);
    Jdark = doDarkJV(eqm_pn.el, voltage_ar, 5);
end
%% 
%Check if the dark eqm solution has a non-zero current 
%i.e., calibrate for errors in numerical integration
 
% voltage_ar = linspace(-0.5, 0.2, 15);
% Jdark = doDarkJV(eqm_PM6.ion, voltage_ar, 5);

%%
sample_point = (length(voltage_ar)-1)/2;
num_points = length(voltage_ar);
dark_current_corr = Jdark.Jvalue - Jdark.Jvalue(voltage_ar == min(abs(voltage_ar)));

figure('Name', 'Dark JV')

semilogy(voltage_ar(Jdark.Jvalue<0), abs(Jdark.Jvalue(Jdark.Jvalue<0)), 'r--')
hold on
semilogy(voltage_ar(Jdark.Jvalue>0), abs(Jdark.Jvalue(Jdark.Jvalue>0)), 'r-')

semilogy(voltage_ar(dark_current_corr<0), abs(dark_current_corr(dark_current_corr<0)), 'b--')
semilogy(voltage_ar(dark_current_corr>=0), abs(dark_current_corr(dark_current_corr>=0)), 'b')

if test_pn == 1
    semilogy(voltage_ar, J_sat_calc*ones(length(voltage_ar)), 'g--')
end

xline(0, 'color', 'black')
yline(0, 'color', 'black')
hold off

xlabel('Voltage (V)')
xlim([voltage_ar(1), voltage_ar(end)])
ylabel('Current Density (A cm^{-2})')
legend({'Uncorrected', '', 'Corrected', '', 'Calculated J_{sat}'})
%ylim([1e-17, 1e-8])