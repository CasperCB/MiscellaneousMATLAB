function [eps] = micaE_IR(wave)
% calculates the infrared dielectric function of muscovite mica based on the optical
% data from Singleton & Shirkley, Appl. Opt., 1983.

% data defined from 400 to 1200 cm^-1

% Clayton B. Casper (11/12/2019)

exp_data = csvread('mica_constants_Singleton_1983.csv'); 

mica_n=interp1(exp_data(:, 1), exp_data(:, 2), wave);
mica_k = interp1(exp_data(:, 1), exp_data(:, 3), wave); 

eps1 = mica_n.^2 - mica_k.^2;
eps2 = 2*mica_n.*mica_k; 

eps = complex(eps1, eps2);

end

