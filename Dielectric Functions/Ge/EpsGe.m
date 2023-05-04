function [eps] = EpsGe(wave, dope)

% function is based off of the Drude model and parameters given in SI of
% Jung, L et. al. ACS Photonics. 2019. 

e = 1.602 * 10^-19; % C
eps0 = 8.854 * 10^-12; % F/m -> (s^2 * C^2)/(m^3 * kg)
c = 2.99792 * 10^8; % m/s
m0 = 9.1095 * 10^-31; % kg

effmass_n = 0.12 * m0; % electron
epsGe = 16; 

freq = 2*pi*c*(wave*100);

wp = sqrt((dope.*(e^2))./(effmass_n*eps0*epsGe));

% gamma = (e./(effmass_n.*mob));
gamma = 116*(10^12); 

eps = epsGe.*(1 - ((wp.^2)./(freq.*(freq + (1i.*gamma))))); 