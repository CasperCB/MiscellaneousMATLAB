function [Alpha_eff] = Alpha_FDM(Trad, L, height, epsSamp, epsSub, angle)
% Finite-dipole model for calculating polarizability, based off 
% of Jung 2019

% added a Fresnel reflection term (10/03/19)
addpath('D:\Atkin Research\Code & Simulation\MATLAB\Tip Dipole Models');
addpath('D:\Atkin Research\Code & Simulation\MATLAB\Dielectric Functions'); 

%% define tip charge distribution
% induced charges
W0 = 1.31*Trad; % position of charge Q0 with respect to radius, m
W1 = 0.5*Trad; % position of charge Q1 with respect to radius, m
g = 0.7*exp(0.06*1i); % contributing portion of near-field induced charge

% finite-dipole geometric factors f0 and f1, formula broken into two parts
f0_1 = (g-((Trad+2*height+W0)./(2*L)));
f0_2 = (log((4*L)./(Trad+(4*height)+(2*W0))))./(log((4*L)./Trad)); 
f0 = f0_1.*f0_2;

f1_1 = (g-((Trad+2*height+W1)./(2*L)));
f1_2 = (log((4*L)./(Trad+(4*height)+(2*W1))))./(log((4*L)./Trad)); 
f1 = f1_1.*f1_2;

%% sample polarizability & reflection
beta = (epsSamp-1)./(epsSamp+1); % sample polarizability term

% substrate refractive index
[n_sub k_sub] = eps2nk(real(epsSub), imag(epsSub)); 

% Fresnel reflection coefficient for p-polarized light
rp = sqrt(1 - (1./n_sub.*sin(angle)).^2)  - n_sub.*cos(angle) ./...
    sqrt(1 - (1./n_sub.*sin(angle)).^2)  + n_sub.*cos(angle);

%% complex polarizability of tip-sample system
Alpha_eff = ((1+rp).^2).*(0.5.*((beta.*f0)./(1-(beta.*f1))))+1; 

end

