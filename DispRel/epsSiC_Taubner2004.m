function [epsSiCTaubner] = epsSiC_Taubner2004(wavenumber, carrConc)

%% constants in SI
eps0 = 8.854 * 10^-12; % F/m -> (s^2 * C^2)/(m^3 * kg)
c = 2.99792 * 10^8; % m/s
e = 1.602 * 10^-19; % C
wave2Hz = 100*c; 

%% parameters defined in paper
% carrConc = 1*(10^15)*(100^3); % m-3
f0 = 788*wave2Hz; % Hz
gamma = 6.8*wave2Hz; % Hz
m_e = 2.88*(10^-31); % kg
tau = 10^-14; % s

%% dielectric function
freq = wavenumber*wave2Hz; % Hz

% sample dielectric function
epsSiCTaubner = 6.49 + (3.23 ./ (1-((freq.^2)./(f0^2))-(1i.*gamma.*freq./(f0.^2)))) - ...
     (carrConc.*(e^2) ./ (eps0.*m_e.*((freq.^2)+(1i.*freq./tau))));
%  
% epsSiCTaubner = 6.49 + 3.23./(1-(freq.^2./f0.^2)-(1i.*gamma.*freq./(f0.^2))) - (carrConc.*(e^2) ./ (eps0.*m_e.*((freq.^2) + 1i.*freq./tau)));
 
end
 