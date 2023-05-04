% This script plots a dispersion relation for a surface plasmon-polartion
% (SPP and/or a surface phonon-polariton (SPhP) via their dielectric functions.

% (Clayton Casper 3/24/17)

% updated on 11/12/19 to be more generalized

addpath(genpath('D:\Atkin Research\Code & Simulation\MATLAB\Tip Dipole Models\')); 
addpath(genpath('D:\Atkin Research\Code & Simulation\MATLAB\Dielectric Functions\')); 

clear; clc;

%% constants
c = 2.998*10^8; % speed of light, m/s
m_e = 9.109*10^-31; % electron mass, kg
e = 1.602*10^-19; % unit of elementary charge, C
eps_0 = 8.85*10^-12; % vacuum permittivity, F/m

%% define parameters
wavelength = linspace(1, 1000000, 1000000).*(10^-9);
wave = (1./(wavelength))./100; 
freq = c./wavelength; 

%% calculate dispersion
eps = epsSiC_Taubner2004(wave, 0); % permittivity of substrate
epsD = 1; % permittivity of surrounding dielectric

% disp = (freq./c).*sqrt((eps.*epsD)./(eps + epsD)); 
disp = (1./wavelength).*sqrt((eps.*epsD)./(eps + epsD)); 

%% plot dispersion
% plot dispersion
hold on
plot((1./wavelength).*1e-6, wave, '-k', 'LineWidth', 1.5);
plot(disp.*1e-6, wave, '-r', 'LineWidth', 1.5); 

% plot Reststrahlen region
% create vectors
% TO_vect = freq_TO.*(ones(6,1));
% LO_vect = freq_LO.*(ones(6,1));
% 
% plot(freq./c, (c./TO_vect).*(10^6), '--k', 'LineWidth', 1.5);
% plot(freq./c, (c./LO_vect).*(10^6), '--k', 'LineWidth', 1.5);

% set(gca, 'XLim', [0.25*(10^5) 2.5*(10^5)]); 
set(gca, 'XLim', [0.084 0.2]); 
set(gca, 'YLim', [850 950]);
set(gca, 'Fontsize', 20); 
xlabel('\bf transverse wavevector \rm (\mum^-^1)'); 
ylabel('\bf wavenumber \rm (cm^-^1)'); 
box on;

% hold off

%% Calculate grating periodicity
% ang = pi/2;
% inc = 10.7*(1e-6); 
% eps1 = epsSiC_Taubner2004((1./(inc))./100, 0);
% out = (1./inc).*sqrt((eps1.*epsD)./(eps1 + epsD)); 
% period = 1/(out-(1/inc)*sin(ang));