% Empirical model for boron-doped Si from Thurber et al. (1980)
% T = 300 K
% Clayton B. Casper (05/14/2021)

function [rho] = Thurber_rho_pSi(dope)

% dope = dopant density, cm^-3

%% fit parameters
qrhoN_min = 0.00215; % (V*sec/cm^2)
qrhoN_max = 0.02053; % (V*sec/cm^2)
Nref = 4.09e18; % cm^-3
alpha = -0.727;

q = 1.602e-19; 
%% qrhoN

qrhoN = qrhoN_min + ((qrhoN_max - qrhoN_min)./ (1 + (dope./Nref).^alpha));

rho = qrhoN./(q*dope); % ohm*cm