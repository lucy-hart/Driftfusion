%The _dark file has a finer point sapcing as I wanted to see if this
%reduced the zero error
%It did not...
% parC60 = pc('Input_files/SAM_MAFACsPbIBr_C60.csv');
% parPM6 = pc('Input_files/SAM_MAFACsPbIBr_PM6Y6.csv');
% parPM7 = pc('Input_files/SAM_MAFACsPbIBr_PM7Y6.csv');
% parPBDBT = pc('Input_files/SAM_MAFACsPbIBr_PBDBTY6.csv');
parC60 = pc('Input_files/SAM_MAFACsPbIBr_C60.csv');
parC60.AbsTol_vsr = 1e-20;
% parY6 = pc('Input_files/SAM_MAFACsPbIBr_Y6.csv');
parPM6 = pc('Input_files/SAM_MAFACsPbIBr_PM6Y6_ShowInterface.csv');
% parPM7 = pc('Input_files/SAM_MAFACsPbIBr_PM7Y6_BHJSurf.csv');
% parPBDBT = pc('Input_files/SAM_MAFACsPbIBr_PCE12Y6_BHJSurf.csv');

eqm_C60 = equilibrate(parC60);
% parPM6.RelTol = 1e-9;
% parPM7.RelTol = 1e-9;
% parPBDBT.RelTol = 1e-9;
%eqm_PM6 = equilibrate(parPM6);
% eqm_PM7 = equilibrate(parPM7);
% eqm_PBDBT = equilibrate(parPBDBT);
% parC60.AbsTol = 1e-12;
% parC60.RelTol = 1e-9;
%%
%Test using pn junction
test_pn = 0;
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
t_hold = 60;

voltage_ar = [-5 -4 -3 -2 -1 -0.5 0 0.1];
% voltage_ar = linspace(-0.5, 0.1, 7);
Jdark = doDarkJV(eqm_C60.el, voltage_ar, t_hold);
%Jdark2 = doDarkJV(eqm_PM6.ion, voltage_ar, t_hold);
% Jdark3 = doDarkJV(eqm_PBDBT.ion, voltage_ar, t_hold);
% Jdark4 = doDarkJV(eqm_PM7.ion, voltage_ar, t_hold);

%%
sample_point = (length(voltage_ar)-1)/2;
num_points = length(voltage_ar);
% dark_current_corr = Jdark.Jvalue - Jdark.Jvalue(voltage_ar == min(abs(voltage_ar)));
% dark_current_corr2 = Jdark2.Jvalue - Jdark2.Jvalue(voltage_ar == min(abs(voltage_ar)));
% dark_current_corr3 = Jdark3.Jvalue - Jdark3.Jvalue(voltage_ar == min(abs(voltage_ar)));
% dark_current_corr4 = Jdark4.Jvalue - Jdark4.Jvalue(voltage_ar == min(abs(voltage_ar)));

semilogplot = 1;

figure('Name', 'Dark JV Raw Data', 'Position', [50 50 1000 1000])

if semilogplot == 1        
        semilogy(voltage_ar(Jdark.Jvalue>0), abs(Jdark.Jvalue(Jdark.Jvalue>0)), 'k-')
        hold on        
        semilogy(voltage_ar(Jdark.Jvalue<0), abs(Jdark.Jvalue(Jdark.Jvalue<0)), 'k--')
%         semilogy(voltage_ar(Jdark4.Jvalue>0), abs(Jdark4.Jvalue(Jdark4.Jvalue>0)), 'Color', [0 0.4470 0.7410])
%         hold on
%         semilogy(voltage_ar(Jdark4.Jvalue<0), abs(Jdark4.Jvalue(Jdark4.Jvalue<0)), 'Color', [0 0.4470 0.7410], 'LineStyle', '--')
        semilogy(voltage_ar(Jdark2.Jvalue>0), abs(Jdark2.Jvalue(Jdark2.Jvalue>0)), 'r-')
        % hold on 
        semilogy(voltage_ar(Jdark2.Jvalue<0), abs(Jdark2.Jvalue(Jdark2.Jvalue<0)), 'r--')        
%         semilogy(voltage_ar(Jdark3.Jvalue>0), abs(Jdark3.Jvalue(Jdark3.Jvalue>0)), 'Color', [0.4660 0.6740 0.1880])
%         semilogy(voltage_ar(Jdark3.Jvalue<0), abs(Jdark3.Jvalue(Jdark3.Jvalue<0)), 'Color', [0.4660 0.6740 0.1880], 'LineStyle', '--')
elseif semilogplot == 0
        hold on
        plot(voltage_ar, Jdark.Jvalue, 'k-')
        plot(voltage_ar, Jdark3.Jvalue, 'Color', [0.4660 0.6740 0.1880])
        plot(voltage_ar, Jdark2.Jvalue, 'r-')
        xline(0, 'color', 'black')

end

yline(0, 'color', 'black')
hold off

set(gca, 'FontSize', 25)
xlabel('Voltage (V)', 'FontSize', 25)
xlim([-5, voltage_ar(end)])
ylim([1e-14, 1])
ylabel('Current Density (A cm^{-2})', 'FontSize', 25)
legend({' CsFAMA', '', ' +0.35 eV', '', ' +0.25 eV', '', ' +0.15 eV', ''}, 'FontSize', 25, 'Location', 'northwest')

% figure('Name', 'Dark JV Corrected Data')
% 
% % semilogy(voltage_ar(dark_current_corr>=0), abs(dark_current_corr(dark_current_corr>=0)), 'k')
% semilogy(voltage_ar(dark_current_corr2<0), abs(dark_current_corr2(dark_current_corr2<0)), 'r--')
% hold on
% % semilogy(voltage_ar(dark_current_corr<0), abs(dark_current_corr(dark_current_corr<0)), 'k--')
% % semilogy(voltage_ar(dark_current_corr2<0), abs(dark_current_corr2(dark_current_corr2<0)), 'r--')
% semilogy(voltage_ar(dark_current_corr2>=0), abs(dark_current_corr2(dark_current_corr2>=0)), 'r-')
% % semilogy(voltage_ar(dark_current_corr3<0), abs(dark_current_corr3(dark_current_corr3<0)), 'Color', [0.4660 0.6740 0.1880], 'LineStyle', '--')
% % semilogy(voltage_ar(dark_current_corr3>=0), abs(dark_current_corr3(dark_current_corr3>=0)), 'Color', [0.4660 0.6740 0.1880])
% 
% % xline(0, 'color', 'black')
% % yline(0, 'color', 'black')
% hold off
% 
% xlabel('Voltage (V)')
% xlim([0, voltage_ar(end)])
% ylabel('Current Density (A cm^{-2})')
% % legend({'C_{60}', 'PM6:Y6'}, 'Location', 'northwest')
