function [epsAuComplex] = GoldE(wavenumber)
%%   Defines Dielectric Constant of Gold from Drude Model
% % 
% e = 1.602 * 10^-19; % C
% eps0 = 8.854 * 10^-12; % F/m -> (s^2 * C^2)/(m^3 * kg)
% c = 2.99792 * 10^8; % m/s
% m0 = 9.1095 * 10^-31; % kg
% cc2cm = 10^-6; % cm^3 -> m^3
% % nAu = (5.9 * 10^22) / cc2cm; % start w/ cm^-3 
% 
% wavelength = (1 ./ wavenumber) .* 10^4; % cm^-1 -> um
% freq = (2 * pi * c) ./ (wavelength .* 10^-6);
% % 
% % wp_Au = sqrt((nAu * e^2)/(m0 * eps0));
% wp_Au = 2043.20100*(10^12); % Olmon, Hz
% gamma_Au = 1./(14*10^15); % Olmon, Hz
% 
% eps_Au = 1 - (wp_Au^2 ./ (freq.^2 + (1i .* freq .* gamma_Au)));
% epsAuComplex = eps_Au;

%% Defines Dielectric Function of Gold from Ellipsometry Data

AuExpData = importdata('Au_dielectric_Evaporated_01.dat', '\t'); % From Raschke Au Simulation
AuWavenumber = 0.01 ./ AuExpData.data(:, 2); % Convert to wavenumber

eps1Au=interp1(AuWavenumber, AuExpData.data(:, 3), wavenumber);
eps2Au=interp1(AuWavenumber, AuExpData.data(:, 4), wavenumber);

epsAuComplex = eps1Au + 1i.*eps2Au;

%% Old Code
% gomegaf=8.06e2;
% gomegat=2.16e2;
% gomegap=6.20e4;
% gEps1=-(gomegap.^2)./(ROIcm.^2+gomegat.^2);
% gEps2=gomegap.^2.*gomegaf./(ROIcm.^3+ROIcm.*gomegat^2);
% 
% epsAuComplex = gEps1 + 1i .* gEps2;

end