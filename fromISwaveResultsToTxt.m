function fromISwaveResultsToTxt(ISwave_results, prefix)
% save the main data from an ISwave_struct created by ISwave_full_exec* to
% txt files, ideally easy to import with Origin (from OriginLab)

%% create header

% check which was the variable being explored
if numel(unique(ISwave_results.Int)) > 1
    legend_text = ISwave_results.Int;
    legend_append = ' sun';
else
    legend_text = ISwave_results.Vdc;
    legend_append = ' Vdc';
end

% round to two significant digits
legendImpedance = round(legend_text, 2, 'significant');
% this will start from 1 sun and go to dark
legendImpedance = string(legendImpedance);
% add sun to numbers in legend
legendImpedance = strcat(legendImpedance, legend_append);
% replace zero in legend with dark
legendImpedance(legendImpedance=="0 sun") = "dark";

headerFrequencyIntVdc = ['Frequency', legendImpedance'];
headerReIm = repelem(strcat(legendImpedance, ' real'), 2);
headerReIm(2:2:end) = legendImpedance;
headerNyquist = ['Frequency', headerReIm'];

%% get measure units

unitsCap = ['Hz', repelem("F/cm\+(2)", length(legendImpedance))];
unitsNyquist = ['Hz', repelem("\g(W)·cm\+(2)", 2*length(legendImpedance))];
unitsZabs = ['Hz', repelem("\g(W)·cm\+(2)", length(legendImpedance))];
unitsPhase = ['Hz', repelem("degrees", length(legendImpedance))];

%% get data

% all the lines are the same for ISwave
frequencies = ISwave_results.Freq(1, :)';

% capacitance
cap = ISwave_results.cap';
dataCap = [frequencies, cap];

% ionic capacitance
capIonic = ISwave_results.cap_idrift';
dataCapIonic = [frequencies, capIonic];

% recombination capacitance
capRec = ISwave_results.cap_U';
dataCapRec = [frequencies, capRec];

% accumulating current capacitance
capAcc = ISwave_results.cap_dQ';
dataCapAcc = [frequencies, capAcc];

% nyquist
impedance_re = ISwave_results.impedance_re';
impedance_im = ISwave_results.impedance_im';
impedance = impedance_im(:, [1;1]*(1:size(impedance_im, 2)));
impedance(:, 1:2:end) = impedance_re;
dataNyquist = [frequencies, impedance];

% absolute impedance
impedance_abs = ISwave_results.impedance_abs';
dataZabs = [frequencies, impedance_abs];

% phase
Zphase = -ISwave_results.J_phase';
dataZphase = [frequencies, Zphase];

% ionic phase
ZphaseIonic = -ISwave_results.J_i_phase';
dataZphaseIonic = [frequencies, ZphaseIonic];

%% join fields

toBeSavedCap = [headerFrequencyIntVdc; unitsCap; dataCap];

toBeSavedCapIonic = [headerFrequencyIntVdc; unitsCap; dataCapIonic];

toBeSavedCapRec = [headerFrequencyIntVdc; unitsCap; dataCapRec];

toBeSavedCapAcc = [headerFrequencyIntVdc; unitsCap; dataCapAcc];

toBeSavedNyquist = [headerNyquist; unitsNyquist; dataNyquist];

toBeSavedZabs = [headerFrequencyIntVdc; unitsZabs; dataZabs];

toBeSavedZphase = [headerFrequencyIntVdc; unitsPhase; dataZphase];

toBeSavedZphaseIonic = [headerFrequencyIntVdc; unitsPhase; dataZphaseIonic];

%% save csv

fid_cap = fopen([prefix '-cap.txt'], 'wt+');
fid_capIonic = fopen([prefix '-cap_ionic.txt'], 'wt+');
fid_capRec = fopen([prefix '-cap_recombination.txt'], 'wt+');
fid_capAcc = fopen([prefix '-cap_accumulating.txt'], 'wt+');
fid_nyquist = fopen([prefix '-nyquist.txt'], 'wt+');
fid_Zabs = fopen([prefix '-Zabs.txt'], 'wt+');
fid_phase = fopen([prefix '-Zphase.txt'], 'wt+');
fid_phaseIonic = fopen([prefix '-Zphase_ionic.txt'], 'wt+');

for i = 1:size(toBeSavedCap, 1)
    fprintf(fid_cap, '%s\t', toBeSavedCap(i, :));
    fprintf(fid_cap, '\n');
    
    fprintf(fid_capIonic, '%s\t', toBeSavedCapIonic(i, :));
    fprintf(fid_capIonic, '\n');
    
    fprintf(fid_capRec, '%s\t', toBeSavedCapRec(i, :));
    fprintf(fid_capRec, '\n');
    
    fprintf(fid_capAcc, '%s\t', toBeSavedCapAcc(i, :));
    fprintf(fid_capAcc, '\n');
    
    fprintf(fid_nyquist, '%s\t', toBeSavedNyquist(i, :));
    fprintf(fid_nyquist, '\n');
    
    fprintf(fid_Zabs, '%s\t', toBeSavedZabs(i, :));
    fprintf(fid_Zabs, '\n');
    
    fprintf(fid_phase, '%s\t', toBeSavedZphase(i, :));
    fprintf(fid_phase, '\n');
    
    fprintf(fid_phaseIonic, '%s\t', toBeSavedZphaseIonic(i, :));
    fprintf(fid_phaseIonic, '\n'); 
end

fclose(fid_cap);
fclose(fid_capIonic);
fclose(fid_capRec);
fclose(fid_capAcc);
fclose(fid_nyquist);
fclose(fid_Zabs);
fclose(fid_phase);
fclose(fid_phaseIonic);