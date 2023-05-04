function [epsPtComplex] = epsPt(wavenumber)
% Defines Dielectric Constant of Platinum

%% constants
% e = 1.602 * 10^-19; % C
% eps0 = 8.854 * 10^-12; % F/m -> (s^2 * C^2)/(m^3 * kg)
c = 2.99792 * 10^8; % m/s
% m0 = 9.1095 * 10^-31; % kg
% cc2cm = 10^-6; % cm^3 -> m^3
% nPt = (5.9 * 10^22) / cc2cm; % start w/ cm^-3 
wave2Hz = 100 * c; 



% %% Defines dielectric function from plasma frequency
% % parameters
% wp_Pt = 4.15 * 10^4; % units 1/(2pi(cm^-1)), from Ordal 1985
% gamma_Pt = 5.58 * 10^2; % units 1/(2pi(cm^-1)), from Ordal 1985 
% 
% wavelength = (1 ./ wavenumber) .* 10^4; % cm^-1 -> um
% freq = c ./ (wavelength .* 10^-6); % Hz
% 
% wp_Pt_Hz = wp_Pt * 2* pi * wave2Hz; % Hz
% gamma_Pt_Hz = gamma_Pt * 2 * pi * wave2Hz; % Hz
% 
% epsPtComplex = 1 - (wp_Pt_Hz^2 ./ (freq.^2 + (1i .* freq .* gamma_Pt_Hz)));

%% Defines dielectric function from literature model
opticsPt = importdata('Pt_Rakic1998.csv', ','); 

wavelength_ref = opticsPt(:, 1); 
wavenumber_ref = 1./(wavelength_ref./10000); 

nPt = interp1(wavenumber_ref, opticsPt(:, 2), wavenumber);
kPt = interp1(wavenumber_ref, opticsPt(:, 3), wavenumber); 

epsPtComplex = (nPt + 1i.*kPt).^2; 
end