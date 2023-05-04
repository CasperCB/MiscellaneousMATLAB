%% this dielectric function is for artifical changes to Si epsilon
%%
% wavenumber = wavenumber in cm^-1
% dope = doping density in m^-3
% field = electric field in V/m
% beta = mobility scalar
% delta = effective mass scalar 

function [eps] = nSiE_toy(wavenumber, dope, field, beta, delta)
e = 1.602 * 10^-19; % C
eps0 = 8.854 * 10^-12; % F/m -> (s^2 * C^2)/(m^3 * kg)
c = 2.99792 * 10^8; % m/s
m0 = 9.1095 * 10^-31; % kg
effmass_n = 0.273 * m0 * delta; % electron, Riffe 2002
epsSi = 11.68; 

freq = 2*pi*c*(wavenumber*100);

wp = sqrt((dope.*(e^2))./(effmass_n*eps0*epsSi));

mob = mobilitySi_Masetti(dope./(10^6), field, true).*beta./(1e4); 
% mob = ((Thurber_mobilitynSi(dope./(10^6)))/(10^4)).*beta; 

gamma = (e./(effmass_n.*mob));

eps = epsSi.*(1 - ((wp.^2)./(freq.*(freq + (1i.*gamma))))); 

end
