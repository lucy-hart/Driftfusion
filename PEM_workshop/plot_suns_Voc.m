function plot_suns_Voc(illuminations, Vocs, output_filename, par)

    log_ill = log(illuminations);
    coeffs = polyfit(log_ill, Vocs, 1);
    m = coeffs(1)/(par.kB*par.T);
    %display(['m = ' num2str(coeffs(1)/(par.kB*par.T))])
    
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

    headers = cell(1, 2);
    headers{1,1} = "Light Intensity";
    headers{1,2} = "Voc";
    values = [illuminations' Vocs'];

    % Create the table and write to file
    output_cell = [headers; num2cell(values)];
    T = cell2table(output_cell);
    writetable(T, ['./Output_files/', output_filename, ' Suns Voc.txt'], 'Delimiter', 'tab', 'WriteVariableNames', 0);

end