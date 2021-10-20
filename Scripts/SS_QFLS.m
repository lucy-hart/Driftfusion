%Script for investigating amount of QFLS at SC and quenching of PL between
%OC and SC as a function of interlayer material and width

%Not inculding ion motion (for now)
%% Start code
tic
initialise_df

%% Read in data and create parameters
par_10 = pc('Input_files/ptaa_mapi_pcbm.csv');

par_50 = par_10;
par_50.d(1) = 5*10^-6;
par_50 = refresh_device(par_50);

par_90 = par_10;
par_90.d(1) = 9*10^-6;
par_90 = refresh_device(par_90);

devices = {par_10, par_50, par_90};

%% Find eqm solutions at 1 sun intensity
eqm_solutions_dark = cell(1,3);
eqm_solutions_light = cell(1,3);
for i = 1:3
    eqm_solutions_dark{i} = equilibrate(devices{i});
    eqm_solutions_light{i} = changeLight(eqm.el,1,1);
end

%% Perform CV scans
CV_solutions = cell(1,3);
for j = 1:3
    sol = eqm_solutions_dark{j};
    CV_solutions{j} = doCV(sol.el, 1, 0, 1.2, 0, 100e-3, 1, 241);
    dfplot.JtotVapp(CV_solutions{j},0)
    hold on
end
hold off
legend('10 nm', '50 nm', '90 nm')
xlim([0, 1.3])

toc
