function [epsSi3N4Complex] = epsSi3N4(wavenumber)
%% Defines dielectric function from literature model
opticsSi3N4 = importdata('Si3N4_Kischkat.csv', ','); 

wavelength_ref = opticsSi3N4(:, 1); 
wavenumber_ref = 1./(wavelength_ref./10000); 

nSi3N4 = interp1(wavenumber_ref, opticsSi3N4(:, 2), wavenumber);
kSi3N4 = interp1(wavenumber_ref, opticsSi3N4(:, 3), wavenumber); 

epsSi3N4Complex = (nSi3N4 + 1i.*kSi3N4).^2; 