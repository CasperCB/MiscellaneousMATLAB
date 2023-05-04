function [epsToy] = epsToy(wave)

%% constants in SI
eps0 = 8.854 * 10^-12; % F/m -> (s^2 * C^2)/(m^3 * kg)
c = 2.99792 * 10^8; % m/s
e = 1.602 * 10^-19; % C
wave2Hz = 100*c; 

omega1 = 548*wave2Hz; % resonance frequency, Hz
omega2 = 795*wave2Hz; % resonance frequency, Hz
eps = 4; % static permittivity
gamma = 11*wave2Hz; % broadening term, Hz

num = eps.*((omega2^2)-(omega1^2)); 

epsToy = eps + num./((omega1.^2)-((wave.*wave2Hz).^2) - (1i.*wave.*gamma.*wave2Hz)); 

end
