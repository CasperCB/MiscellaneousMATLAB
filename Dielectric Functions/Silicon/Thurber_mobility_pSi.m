% Empirical model for boron-doped Si from Thurber et al. (1980)
% T = 300 K
% Clayton B. Casper (05/14/2021)

function [mob] = Thurber_mobility_pSi(rho)

% rho = resistivity, Ohm*cm

%% fit parameters
A = 51.6; % cm^2/(V*s)
rho_c = 0.004063; % Ohm*cm
mu_max = 467.3; % cm^2/(V*s)
rho_ref = 0.0794; % Ohm*cm
alpha = -0.808;

%% mu

mob = A.*exp(-rho/rho_c) + mu_max./(1 + (rho./rho_ref).^alpha); 