function [eps] = SiO2E_vis(wavelength)
% calculates the dielectric function of gold based on the optical
% data from Rodriguez-de Marcos

% data defined from 300 to 1510 nm

% input wavelength in nm

% Clayton B. Casper (3/3/2020)

exp_data = csvread('Rodriguez-deMarcos.csv'); 

n = interp1(exp_data(:, 1), exp_data(:, 2), wavelength.*(10^-9));
k = interp1(exp_data(:, 1), exp_data(:, 3), wavelength.*(10^-9)); 

eps1 = n.^2 - k.^2;
eps2 = 2*n.*k; 

eps = complex(eps1, eps2);

end

