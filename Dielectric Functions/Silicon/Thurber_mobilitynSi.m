% Empirical model for phosphorus-doped Si from Thurber et al. (1980)
% T = 300 K
% Clayton B. Casper (10/14/19)

function [mobility] = Thurber_mobilitynSi(dope)

% dope = electrically-active dopant density, cm^-3

%% fit parameters
mu0 = 1; % cm^2/(V*s)
n0 = 10^16; % ohm*cm

A0 = 3.0629;
A1 = -2.2522;
A2 = 0.62327;
A3 = -0.060415;
B1 = -0.69851;
B2 = 0.19716; 
B3 = -0.019950; 

%% mobility

X = log10(dope/n0); 

mulog = (A0 + A1.*X + A2.*(X.^2) + A3.*(X.^3)) ./ (1 + B1.*X + B2.*(X.^2) + B3.*(X.^3)); 

mobility = mu0.*(10.^(mulog)); 