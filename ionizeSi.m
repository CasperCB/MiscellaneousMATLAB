function [N_ion] = ionizeSi(R, N, type)
% This function is meant to calculate the amount of ionized dopants
% in a Si:NW given an encoded dopant concentration. This function 
% is based off of a discussion in Schmidt, V. et al. Adv. Mater. 2009.

% Clayton B. Casper (10/09/19)

% R = nanowire radius in nm
% N = dopant concentration in cm^-3
% type = 1 for N = Nd or 0 for N = Na

%% constants
k = 1.38*10^-23; % J/K, Boltzmann's constant
e = 1.602*10^-19; % C, elementary charge
m = 9.109*10^-31; % kg, mass of an electron
h = 6.626*10^-34; % J*s; Planck's constant

%% parameters
T = 300; % K, temperature

eps_s = 11.67; % static permittivity of Si
Ei0 = 45; % meV, bulk ionization energy of dopants (Diarra 2007)
m_e = 0.273*m; % kg, electron effective mass
m_h = 0.367*m; % kg, hole effective mass

eps_o = 1; % permittivity of surrounding medium, assuming air

%% calculate ionized dopant concentration
Nc = ((2*((2*pi*m_e*k*T)/(h^2)))^(3/2))./(10^6); % cm^-3, conduction band effective DOS
Nv = ((2*((2*pi*m_h*k*T)/(h^2)))^(3/2))./(10^6); % cm^-3, valence band effective DOS

x = (eps_s/eps_o); % permittivity ratio

% empirical function in Ei radius dependence, units eV*nm (from Niquet 2006)
F = (200.674 + 175.739*x + 17.395*(x^2) + 0.0949*(x^3)) / (219.091 + 50.841*x + (x^2)); 

% radius dependent ionization energy from Diarra 2007
Ei = (Ei0/10^3) + (2./(R.*eps_s)).*((eps_s - eps_o)/(eps_s + eps_o)).*F; 

% ionized dopant concentration

if type == 1
    N_ion = (Nc/4).*(exp(-Ei*e/(k*T))).*(-1 + sqrt(1 + (8*(N./Nc).* (exp(Ei*e/(k*T)))))); 
    
else 
    N_ion = (Nv/8).*(exp(-Ei*e/(k*T))).*(-1 + sqrt(1 + (16*(N./Nv).* (exp(Ei*e/(k*T)))))); 
end

end

