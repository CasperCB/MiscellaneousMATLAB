function [eps] = GoldE_vis(wavelength)
% calculates the dielectric function of gold based on the optical
% data from Johnson & Christy

% data defined from 188 to 1940 nm

% input wavelength in nm

% Clayton B. Casper (3/3/2020)

exp_data = csvread('OpticalDataAuJohnson&Christy1972.csv'); 

n = interp1(exp_data(:, 2), exp_data(:, 3), wavelength.*(10^-9));
k = interp1(exp_data(:, 2), exp_data(:, 4), wavelength.*(10^-9)); 

eps1 = n.^2 - k.^2;
eps2 = 2*n.*k; 

eps = complex(eps1, eps2);

end

