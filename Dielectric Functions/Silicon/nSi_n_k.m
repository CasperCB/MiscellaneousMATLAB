% Quick script for calculating the refractive indices of n-Si at a given
% wavelength.

% Clayton B. Casper (05/20/2021)

clear; clc;

%% constants
e = 1.602 * 10^-19; % C, elementary charge constant
eps0 = 8.854 * 10^-12; % F/m -> (s^2 * C^2)/(m^3 * kg), vacuum permittivity
c = 2.99792 * 10^8; % m/s, speed of light
m0 = 9.1095 * 10^-31; % kg, mass of an electron

%% Silicon parameters
effmass_n = 0.273 * m0; % kg, electron effective mass, Riffe 2002
epsSi = 11.68; % dielectric constant

dope = logspace(17, 20, 1000); % cm^-3, dopant density 
wave = 940; % cm^-1, wavenumber 

%% conversions
freq = 2*pi*c*(wave*100); % Hz, wavenumber -> frequency
dope_m3 = dope.*1e6; % cm^-3 -> m^-3

%% permittivity calculation (Drude model)
wp = sqrt((dope_m3.*(e^2))./(effmass_n*eps0*epsSi)); % Hz, plasma frequency

mob = mobilitySi_Masetti(dope, 0, true)/(1e4); % cm^2/(V*s), electron mobility from Masetti model

gamma = (e./(effmass_n.*mob)); % Hz, damping rate

eps = epsSi.*(1 - ((wp.^2)./(freq.*(freq + (1i.*gamma))))); % complex permittivity 

n = sqrt((sqrt((real(eps).^2) + (imag(eps).^2)) + real(eps))./2); % refractive index
k = sqrt((sqrt((real(eps).^2) + (imag(eps).^2)) - real(eps))./2); % extinction coefficient


