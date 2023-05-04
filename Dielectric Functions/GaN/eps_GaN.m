function [eps] = epsGaN(wave, n, dope)
% wave = wavenumber in cm^-1
% n = carrier density in cm^-3
% dope = dopant type, 1 = n-type, 0 = p-type

%% constants
c = 2.998*10^8; % speed of light, m/s
m_e = 9.109*10^-31; % electron mass, kg
e = 1.602*10^-19; % unit of elementary charge, C
eps_0 = 8.85*10^-12; % vacuum permittivity, F/m

%% material parameters
eps_inf = 

%% calculation

wavelength = (1./(wave.*100)); % m, wavelength 
freq = c*100*wave; % Hz, frequency

wp = sqrt((dope*10^6)*(e^2)/(m_eff*eps_0)); 
gamma_e = e/(mu*m_eff); 

eps_e = eps_inf.*(1-((wp^2)./((freq.^2) + 1i.*freq.*gamma_e))); % contribution from free carriers
eps_e2 = -eps_inf.*((wp^2)./((freq.^2) + 1i.*freq.*gamma_e));

eps_ph = eps_inf.*(1 + ((w_LO^2)-(w_TO^2))./((w_TO^2)-(freq.^2)-1i.*freq.*gamma_ph)); % contribution from phonons

eps = eps_e2 + eps_ph; % GaN dielectric function
end

