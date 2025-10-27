function plot_suns_Voc(illuminations, Vocs, par)

    log_ill = log(illuminations);
    coeffs = polyfit(log_ill, Vocs, 1);
    display(['m = ' num2str(coeffs(1)/(par.kB*par.T))])
    
    figure('Name', 'SunsVoc')
    
    semilogx(illuminations, Vocs, 'Marker', 'o', 'LineStyle', 'none')
    hold on
    semilogx(illuminations, coeffs(2) + coeffs(1)*log_ill)
    hold off 
    
    xlim([illuminations(1), illuminations(end)])
    ylim([0.98*Vocs(1), 1.01*Vocs(end)])
    xlabel('Illumination (Suns)')
    ylabel('V_{OC}')
    legend('Data', 'Fit', 'Location', 'southeast')

end