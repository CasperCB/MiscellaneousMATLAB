function [epsSiCAmarie] = epsSiC_OcelicDiss(wavenumber, pol)

%% constants in SI
eps0 = 8.854 * 10^-12; % F/m -> (s^2 * C^2)/(m^3 * kg)
c = 2.99792 * 10^8; % m/s
e = 1.602 * 10^-19; % C
wave2Hz = 100*c; 

%% parameters defined in paper
if pol == 0 % perpendicular to c-axis
    wLO = 971; % cm^-1
    wTO = 797; % cm^-1
    eps = 6.56; 
    gamma_ph = 6.6; % cm^-1
    gamma_e = 450; % cm^-1
    wp = 275; % cm^-1
end

if pol == 1 % parallel to c-axis
    wLO = 967; % cm^-1
    wTO = 782; % cm^-1
    eps = 6.56; 
    gamma_ph = 6.6; % cm^-1
    gamma_e = 450; % cm^-1
    wp = 220; % cm^-1
end

%% dielectric function
% phonon term
eps_ph = ((wLO^2)-(wTO^2))./((wTO^2)-(wavenumber.^2)-(1i*wavenumber*gamma_ph)); 

% plasma term
eps_p = (wp^2)./(-(wavenumber.^2)-(1i*wavenumber*gamma_e)); 

epsSiCAmarie = eps.*(1 + eps_ph + eps_p); 

end
 