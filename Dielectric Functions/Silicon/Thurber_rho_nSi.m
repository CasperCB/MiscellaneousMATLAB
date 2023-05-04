% Empirical model for phosphorus-doped Si from Thurber et al. (1980)
% T = 300 K
% Clayton B. Casper (10/14/19)

function [rho] = Thurber_rho_nSi(N)

% dope = dopant density, cm^-3

%% fit parameters
mu0 = 1; % cm^2/(V*s)
N0 = 10^16; % cm^-3
rho0 = 1; % ohm*cm

A0 = -3.0652;
A1 = 2.1853;
A2 = -0.61080;
A3 = 0.056189;
B1 = -0.67642;
B2 = 0.19542; 
B3 = -0.018100; 

%% mobility

X = log10(N/N0); 

rholog = (A0 + A1.*X + A2.*(X.^2) + A3.*(X.^3)) ./ (1 + B1.*X + B2.*(X.^2) + B3.*(X.^3)); 

rho = (rho0.*(10.^(rholog)))./((1.602*10^-19)*N); 