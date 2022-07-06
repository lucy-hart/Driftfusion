% Function LoadRefrIndex
% This function returns the complex index of refraction spectra, ntotal, for the
% material called 'name' for each wavelength value in the wavelength vector
% 'wavelengths'.  The material must be present in the index of refraction
% library 'Index_of_Refraction_library.xls'.  The program uses linear
% interpolation/extrapolation to determine the index of refraction for
% wavelengths not listed in the library.
function [n_interp, k_interp] = LoadRefrIndex(name,wavelengths)

%Data in IndRefr
IndRefr=readtable('Index_of_Refraction_library.xls', 'VariableNamingRule', 'preserve');

% Load index of refraction data in spread sheet, will crash if misspelled
file_wavelengths=IndRefr.(1);
n=IndRefr.(strcat(name,'_n'));
k=IndRefr.(strcat(name,'_k'));

% Interpolate/Extrapolate data linearly to desired wavelengths
n_interp=interp1(file_wavelengths, n, wavelengths, 'linear', 'extrap');
k_interp=interp1(file_wavelengths, k, wavelengths, 'linear', 'extrap');

%Return interpolated complex index of refraction data
%ntotal = n_interp+1i*k_interp; 

end