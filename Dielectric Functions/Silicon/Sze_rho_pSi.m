function [rho] = Sze_rho_pSi(dope)
% calculates resistivity of p-type Si based on figure in Sze 1986, which
% was digitized using webplotdigitizer

% Clayton B. Casper (12/10/19)

%% import data
data = csvread('Sze_1986_ptype_rho.csv'); 

n_ref = data(:,1); % cm^-3
rho_ref = data(:,2); % ohm*cm

%% interpolation
rho = interp1(n_ref, rho_ref, dope); 

end

