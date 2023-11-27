%Use this file to symmetrically sweep HOMO offset vs k2

%TURN SAVE OFF TO START OFF WITH (final cell)
%% Define parameter space
Donor_HOMO = linspace(-5.3, -5.0, 4);
n_HOMOs = length(Donor_HOMO);
t_hold = 60;
voltage_ar = [-5 -0.5 0];

%%
error_log = zeros(n_HOMOs,1);
soleq = cell(n_HOMOs,1);
solCV = cell(n_HOMOs,1);
results = cell(n_HOMOs,1);
Jdark = cell(n_HOMOs,1);

%% Do (many) JV sweeps
%Remeber to update the work function values if you change these parameters
%between files 
illumination = 1;

%Set this to one if you are using the '..._ShowInterface.csv' file
surface = 0;

%Reset the electrode work functions in each loop to be safe as they are
%changed for the cases where E_LUMO (E_HOMO) is far below (above) the CB
%(VB)
for i = 1:n_HOMOs
    disp(["Donor_HOMO = ", num2str(Donor_HOMO(i)), " eV"])
    
    %par = pc('Input_files/SAM_MAFACsPbIBr_PM6Y6_BHJ_SRH.csv');
    par = pc('Input_files/SAM_MAFACsPbIBr_PM6Y6.csv');
%     par = pc('Input_files/SAM_MAFACsPbIBr_PM6Y6_ShowInterface.csv');

    %Donor HOMO Energetics
    par.Phi_IP(4:6) = Donor_HOMO(i);
    par.EF0(4:6) = (par.Phi_IP(5)+par.Phi_EA(5))/2;
    par.Et(4:6) = (par.Phi_IP(5)+par.Phi_EA(5))/2;
    if surface == 1
        par.Et2(4:6) = par.Phi_IP(5) + 0.05; 
    end

    %par.RelTol = 1e-9;

    par = refresh_device(par);

    soleq{i} = equilibrate(par);

    Jdark{i} = doDarkJV(soleq{i}.ion, voltage_ar, t_hold);
    
    Voc_max = 1.2;
    num_points = 281; 
    while Voc_max >= 1.05
        try
            solCV{i} = doCV(soleq{i}.ion, illumination, -0.2, Voc_max, -0.2, 25e-3, 1, num_points);
            error_log(i) = 0;
            results{i} = CVstats(solCV{i});
            Voc_max = 0;
        catch
            if Voc_max > 1.05
                warning("Ionic JV solution failed, reducing Vmax by 0.03 V")
                Voc_max = Voc_max - 0.03;
                num_points = num_points - 6;
            elseif Voc_max == 1.05
                warning("Ionic JV solution failed.")
                error_log(i) = 1;
                results{i} = 0;
            end
        end
    end

end


%% Plot results 
Stats_array = zeros(n_HOMOs, 4);
for i = 1:n_HOMOs
    try
        Stats_array(i,1) = 1e3*results{i}.Jsc_f;
        Stats_array(i,2) = results{i}.Voc_f;
        Stats_array(i,3) = results{i}.FF_f;
        Stats_array(i,4) = results{i}.efficiency_f;
    catch
        warning('No Stats')
        Stats_array(i,:) = 0;
    end
end


%% Run JV for C60 device to get the Jsc value 
parC60 = pc('Input_files/SAM_MAFACsPbIBr_C60.csv');
eqm_QJV_C60 = equilibrate(parC60);
CV_sol_C60 = doCV(eqm_QJV_C60.ion, illumination, -0.2, 1.2, -0.2, 25e-3, 1, 281);
stats_C60 = CVstats(CV_sol_C60);

%%
figure('Name', 'Jsc and Jd vs Injection Barrier', 'Position', [50 50 1300 1500])
ax = axes;
yyaxis('left')
yyaxis('right')
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = [1 0 0];
ax.LineWidth = 2;

box on 
hold on
yyaxis left
set(gca, 'Fontsize', 25)
ylabel('J_{SC} (mA cm^{-2})', 'FontSize', 30)
ylim([14.5, 25])
plot(Donor_HOMO + 5.5, -Stats_array(:, 1), 'color', 'black', 'LineStyle', 'none', 'Marker', 'square', 'MarkerSize', 15, ...
    'MarkerFaceColor', 'black')
yline(abs(1e3*stats_C60.Jsc_f), 'Color', 'black', 'LineStyle', '--', 'LineWidth', 3)
txt = 'J_{SC,ref}';
text(0.21, 21.2, txt, 'FontSize', 30)

yyaxis right
set(gca, 'Fontsize', 25, 'YScale', 'log')
ylabel('J_{d} (A cm^{-2})', 'FontSize', 30)
ylim([9e-10, 1e-4])
for i = 1:n_HOMOs
    plot(Donor_HOMO(i)+5.5, abs(Jdark{i}.Jvalue(1)), 'color', 'red', 'LineStyle', 'none', 'Marker', 'square', ...
        'MarkerSize', 15, 'MarkerFaceColor', 'red')
    errorbar(Donor_HOMO(i)+5.5, abs(Jdark{i}.Jvalue(1)), abs(Jdark{i}.Jvalue(2)), 'color', 'red')
end
hold off

xlabel('Energetic Barrier to Hole Injection (eV)', 'FontSize', 30)
xticks([0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5])
xticklabels({'0.20', '0.25', '0.30', '0.35', '0.40', '0.45', '0.50'})
xlim([0.2, 0.5])

ax1 = gca;
save_fig = 0;
if save_fig == 1
    exportgraphics(ax1, ...
    'C:\Users\ljh3218\OneDrive - Imperial College London\PhD\Davide_OrganicPeroHybrid\FigSxa-Minus5VJd.png', ...
    'Resolution', 300)
end

%% Save results and solutions
save_file = 0;
if save_file == 1
    filename = 'PeroBHJ_QuasiTandem_DonorHOMO_Jd_Jsc.mat';
    save(filename, 'solCV')
end
